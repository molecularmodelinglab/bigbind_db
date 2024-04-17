# imports :)
import subprocess
import time
import sys
import os
from glob import glob
import sqlite3
import pandas as pd
from tqdm import tqdm
import urllib.request
from pyfaidx import Fasta
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem import rdDistGeom
import rdkit
from rdkit.Chem import AllChem
from pathlib import Path
from bigbind.db import create_connection
from rdkit.Chem import Descriptors
import yaml
import multiprocessing as mp
from functools import partial
import logging
import signal 
from concurrent import futures
from bigbind.config import CONFIG


log = logging.getLogger(__name__)

#classes-------------------


class timeout:
    def __init__(self, seconds, error_message="Timeout"):
        self.seconds = seconds
        self.error_message = error_message

    def handle_timeout(self, signum, frame):
        raise TimeoutError(self.error_message)

    def __enter__(self):
        signal.signal(signal.SIGALRM, self.handle_timeout)
        signal.alarm(self.seconds)

    def __exit__(self, type, value, traceback):
        signal.alarm(0)


def untar_chembl(out_path, chembl_filename):
    print(f"Untarring {chembl_filename}")
    out_dir = os.path.dirname(out_path)
    cmd = f"tar -xf {chembl_filename} --one-top-level={out_dir}"
    print(f"Running {cmd}")
    subprocess.run(cmd, shell=True, check=True)
    # want to make this work for future chembl versions as well
    db_file = glob(os.path.join(out_dir, "chembl_*/chembl_*_sqlite/chembl_*.db"))[0]
    os.rename(db_file, out_path)
    return out_path


def get_chembl_con(chembl_db_file):
    """Gets the connection to the chembl sqlite database"""
    con = sqlite3.connect(chembl_db_file)
    return con


def get_crossdocked_chembl_activities(csv_path, con, prev_output=None):
    """Get all activities (with some filtering for quality) of small
    molecules binding to proteins whose structures are in the crossdocked
    dataset"""

    if prev_output is not None:
        return prev_output

    query = f"""

    SELECT md.chembl_id AS compound_chembl_id,
    cs.canonical_smiles,
    act.standard_type,
    act.standard_relation,
    act.standard_value,
    act.standard_units,
    act.pchembl_value,
    act.potential_duplicate,
    COALESCE(act.data_validity_comment, 'valid') as data_validity_comment,
    a.confidence_score,
    td.chembl_id AS target_chembl_id,
    td.target_type,
    c.accession as protein_accession,
    a.chembl_id as assay_id
    FROM target_dictionary td
    JOIN assays a ON td.tid = a.tid
    JOIN activities act ON a.assay_id = act.assay_id
    JOIN molecule_dictionary md ON act.molregno = md.molregno
    JOIN compound_structures cs ON md.molregno = cs.molregno
    JOIN target_type tt ON td.target_type = tt.target_type
    JOIN target_components tc ON td.tid = tc.tid
    JOIN component_sequences c ON tc.component_id = c.component_id
    AND tt.parent_type  = 'PROTEIN' 
    AND act.pchembl_value IS NOT NULL
    """
    with open(csv_path, "w"):
        pass

    approx_tot = 3212150
    chunksize = 1000
    chunks = pd.read_sql_query(query, con, chunksize=chunksize)
    header = True
    for i, chunk in enumerate(tqdm(chunks, total=int(approx_tot / chunksize))):
        chunk.to_csv(csv_path, header=header, mode="a", index=False)

        header = False

    return csv_path


def reporthook(count, block_size, total_size):
    global start_time
    if count == 0:
        start_time = time.time()
        return
    duration = time.time() - start_time
    progress_size = int(count * block_size)
    speed = int(progress_size / (1024 * duration))
    percent = int(count * block_size * 100 / total_size)
    sys.stdout.write("\r...%d%%, %d MB, %d KB/s, %d seconds passed" %
                    (percent, progress_size / (1024 * 1024), speed, duration))
    sys.stdout.flush()

def download_chembl(desired_db_path, desired_csv_path):
    # Check if the CSV file exists
    if not os.path.isfile(desired_csv_path):
        print(f'{desired_csv_path} does not exist, downloading.')
        chembl_tarred = urllib.request.urlretrieve(
            f"https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_{CONFIG.chembl_version}_sqlite.tar.gz",
            "./temp_data/chembl_33_sqlite.tar.gz",
            reporthook
        )[0]
        print(f'{desired_csv_path} downloaded')
        chembl_db_file = untar_chembl(desired_db_path, chembl_tarred)

        con = get_chembl_con(chembl_db_file)

        chembl_csv_name = get_crossdocked_chembl_activities(desired_csv_path, con)

    # Read the CSV file into a DataFrame
    dataframe = pd.read_csv(desired_csv_path)

    return dataframe


# takes in chembl dataframe and outputs molecules

def gotzinc(molecule):
    yayatoms = ["H", "C", "N", "O", "F", "S", "P", "Cl", "Br", "I"]
    atoms = molecule.GetAtoms()
    for i in atoms:
        if i.GetSymbol() not in yayatoms:
            return False
    return True


def get_conformer(mol, chemblid):
    try:
        #times out after 30 seconds when creating conformers
        with timeout(30):
            try:
                p = Path(f"data/molecules/conformers/{chemblid}.sdf")
                if p.is_file():
                    return True

                mol = Chem.AddHs(mol)
                conformers = AllChem.EmbedMultipleConfs(mol, numConfs=1)

                writer = Chem.SDWriter(str(p.resolve()))
                for cid in range(mol.GetNumConformers()):
                    writer.write(mol, confId=cid)
                return True
            except RuntimeError:
                print("unable to work")
                return False
    except TimeoutError:
        print("timedout")
        return False


def molecules_sequence_chunk(molecules):
    for index, row in molecules.iterrows():
        cur_mol = Chem.MolFromSmiles(row["canonical_smiles"])

        # get molecular weight
        molecules.at[index, "molecular_weight"] = Descriptors.ExactMolWt(cur_mol)
        # does it ONLY have elements in the list yayatoms
        molecules.at[index, "zinc_elements"] = gotzinc(cur_mol)

        # create conformer
        #if not os.path.exists("data/molecules/conformers"):
        #    os.makedirs("data/molecules/conformers")
        Path("data/molecules/conformers").mkdir(parents=True, exist_ok=True)

        molecules.at[index, "has_conformer"] = get_conformer(
            cur_mol, row["compound_chembl_id"]
        )
    return molecules

#@task
def create_molecules(chembl_df, break_num):
    molecules = pd.DataFrame()
    molecules = molecules.assign(compound_chembl_id=chembl_df["compound_chembl_id"])
    molecules = molecules.assign(canonical_smiles=chembl_df["canonical_smiles"])

    # make empty columns
    molecules["molecular_weight"] = -1.0
    molecules["zinc_elements"] = False
    molecules["has_conformer"] = False

    #making it the size of break num
    molecules = molecules[:break_num]

    # for debugging
    shape = chembl_df.shape[0]

    if not os.path.exists(os.path.dirname("data/molecules/confomers")):
        os.makedirs(os.path.dirname("data/molecules/confomers"))

    n_jobs = mp.cpu_count() // 2  # Poolsize
    pool = mp.Pool(n_jobs)
    print('Starting Molecule MP')
    chunksize = len(molecules)//n_jobs
    chunks = [molecules[i:i+chunksize] for i in range(0,len(molecules),chunksize)]
    print("loading result")
    result = pool.map(molecules_sequence_chunk, chunks)
    molecules = pd.concat(result)
    print("closing pool")
    pool.close()
    print("joining pool")
    pool.join()

    return molecules


def download_uniprot_fasta(out_path, downloadurl):
    if os.path.exists(out_path.rsplit(".", 1)[0]):  # file already exists
        return out_path
    if not os.path.exists(os.path.dirname(out_path)):
        os.makedirs(os.path.dirname(out_path))
    urllib.request.urlretrieve(downloadurl, out_path, reporthook)
    cmd = f"gzip -d {out_path}"
    print(f"Running {cmd}")
    subprocess.run(cmd, shell=True, check=True)
    return out_path


def extract_id(header):
    return header.split("|")[1]

def proteins_sequence_chunk(proteins, sequences_complete, sequences_uncomplete):
    for index, row in proteins.iterrows():

        if row["uniprotID"] in sequences_complete:
            proteins.at[index, "protein_sequence"] = str(sequences_complete[row["uniprotID"]])

        else:
            proteins.at[index, "protein_sequence"] = str(sequences_uncomplete[row["uniprotID"]])

    return proteins

#@task
def create_proteins(chembl_df, break_num):
    # download uniprot sequences
    print("Creating proteins...")
    download_uniprot_fasta(
        "data/uniprot/uniprot_sprot.fasta.gz",
        "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz",
    )
    print("finished reviewed uniprot")
    download_uniprot_fasta(
        "data/uniprot/uniprot_trembl.fasta.gz",
        "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz",
    )
    print("finished unreviewed uniprot")
    # make proteins Dataframe
    proteins = pd.DataFrame()
    # get ID's
    proteins = proteins.assign(uniprotID=chembl_df["protein_accession"])
    # remove duplicates
    proteins = proteins.drop_duplicates()
    #make it as big as our break_num
    proteins = proteins[:break_num]

    # add new column
    proteins["protein_sequence"] = ""
    print("loading fasta...")
    sequences_complete = Fasta(
        "data/uniprot/uniprot_sprot.fasta", key_function=extract_id
    )
    sequences_uncomplete = Fasta(
        "data/uniprot/uniprot_trembl.fasta", key_function=extract_id
    )

    # Job parameters
    n_jobs = mp.cpu_count() // 2  # Poolsize
    #pool = mp.Pool(n_jobs)
    print('Starting Protein MP')
    chunksize = len(proteins)//n_jobs

    chunks = [proteins[i:i+chunksize] for i in range(0,len(proteins),chunksize)]

    process_partial = partial(proteins_sequence_chunk, sequences_complete=sequences_complete, sequences_uncomplete=sequences_uncomplete)
    with futures.ThreadPoolExecutor(n_jobs) as executor:
        result = executor.map(process_partial, chunks)
    
    
    proteins = pd.concat(result)


    return proteins


#@task
def create_activites(chembl_df, break_num):
    activities = pd.DataFrame()
    activities = activities.assign(
        # molecules for referance
        compound_chembl_id=chembl_df["compound_chembl_id"],
        canonical_smiles=chembl_df["canonical_smiles"],
        # uniprotid protein reference
        proteinID=chembl_df["protein_accession"],
        # all the info related to molecule/protein
        standard_type=chembl_df["standard_type"],
        standard_relation=chembl_df["standard_relation"],
        standard_value=chembl_df["standard_value"],
        standard_units=chembl_df["standard_units"],
        pchembl_value=chembl_df["pchembl_value"],
    )

    activities = activities[:break_num]
    return activities


#@flow
def load_chembl():
    print("Starting Main")

    max_table_len = CONFIG.max_table_len

    df = download_chembl("data/chembl/chembl.db", "data/chembl/chembl.csv")
    print("Loading molecules...")
    molecules = create_molecules(df, max_table_len)
    print("Saving molecules")
    # molecules.to_csv("molecules.csv", index=False)

    print("Loading proteins...")
    proteins = create_proteins(df, max_table_len)
    # proteins.to_csv("proteins.csv", index=False)

    print("Loading activities...")
    activites = create_activites(df, max_table_len)
    # activites.to_csv("activites.csv", index=False)

    con = create_connection()
    molecules.to_sql(con=con, name='molecules', schema='SCHEMA', index=False, if_exists='append')
    proteins.to_sql(con=con, name='proteins', schema='SCHEMA', index=False, if_exists='append')
    activites.to_sql(con=con, name='activites', schema='SCHEMA', index=False, if_exists='append')
    print("done")

if __name__ == "__main__":
    load_chembl()
