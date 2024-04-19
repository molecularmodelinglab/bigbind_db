# Example on how loading data into BigBind should work
# Basically put functions into prefect tasks, especially if
# they can be called concurrently

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
        cur_mol = Chem.MolFromSmiles(row["smiles"])

        # get molecular weight
        molecules.at[index, "molecular_weight"] = Descriptors.ExactMolWt(cur_mol)
        # does it ONLY have elements in the list yayatoms
        molecules.at[index, "zinc_elements"] = gotzinc(cur_mol)

        # create conformer
        #if not os.path.exists("data/molecules/conformers"):
        #    os.makedirs("data/molecules/conformers")
        Path("data/molecules/conformers").mkdir(parents=True, exist_ok=True)

        molecules.at[index, "has_conformer"] = get_conformer(
            cur_mol, row["chembl_id"]
        )
    return molecules

#@task
def create_molecules(chembl_df, break_num):
    molecules = pd.DataFrame()
    molecules = molecules.assign(chembl_id=chembl_df["compound_chembl_id"])
    molecules = molecules.assign(smiles=chembl_df["canonical_smiles"])
    #dropping duplicates 
    molecules = molecules.drop_duplicates()

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
    
    molecules["chembl_id"] = [ int(''.join(c for c in x if c.isdigit())) for x in molecules["chembl_id"]]
    

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

        if row["uniprot"] in sequences_complete:
            proteins.at[index, "sequence"] = str(sequences_complete[row["uniprot"]])

        elif sequences_uncomplete is not None and row["uniprot"] in sequences_uncomplete:
            proteins.at[index, "sequence"] = str(sequences_uncomplete[row["uniprot"]])
        else:
            proteins.at[index, "sequence"] = "none"

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
    proteins = proteins.assign(uniprot=chembl_df["protein_accession"])
    # remove duplicates
    proteins = proteins.drop_duplicates()
    #make it as big as our break_num
    proteins = proteins[:break_num]

    # add new column
    proteins["sequence"] = ""
    print("loading fasta...")
    sequences_complete = Fasta(
        "data/uniprot/uniprot_sprot.fasta", key_function=extract_id
    )
    if not CONFIG.only_use_complete_uniprots:
        sequences_uncomplete = Fasta(
            "data/uniprot/uniprot_trembl.fasta", key_function=extract_id
        )
    else:
        sequences_uncomplete = None

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
    
    #making sure its in form for sql table
    # proteins["id"] = [ int(''.join(c for c in x if c.isdigit())) for x in proteins["id"]]
    # proteins = proteins.drop_duplicates(subset=["id"])
    proteins = proteins.drop_duplicates(subset=["sequence"])

@task
def download_chembl():
    """ Download the ChEMBL database into CONFIG.temp_data
    and return the filename of the SQLite file.
    (I imagine prefect will have issues passing around
    SQLite connections lol) """
    raise NotImplementedError

@task
def add_chembl_molecules(chembl_filename):
    raise NotImplementedError

@task
def add_chembl_proteins(chembl_filename):
    raise NotImplementedError

@task
def add_chembl_activities(chembl_filename):
    raise NotImplementedError

@flow
def load_chembl():
    chembl_filename = download_chembl()
    mol_future = add_chembl_molecules.submit(chembl_filename)
    prot_future = add_chembl_proteins.submit(chembl_filename)
    # wait for both to finish
    mol_future.result()
    prot_future.result()
    # now we can add the activities
    add_chembl_activities(chembl_filename)