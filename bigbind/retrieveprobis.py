#testing multiple pdb, cen files
#Retrieve PDB ID, Chain ID, Residue ID for residues <x angstroms from binding site centroid for each chain in protein

import numpy as np
import pandas as pd
import warnings
import Bio
from Bio import Align
from Bio.Seq import Seq
from Bio.PDB import *
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from glob import glob
from bigbind.config import CONFIG
from bigbind.probis_tables import add_example_table, add_data_to_example
from bigbind.db import create_connection
parser = PDBParser(PERMISSIVE=True)
warnings.simplefilter('ignore', PDBConstructionWarning)
aligner = Align.PairwiseAligner()

#probisdf = pd.DataFrame(columns=['PDB ID', 'Chain ID', 'Binding site residues'])
#code for sql table creation
con = create_connection()
probis_tbl = add_example_table(con)
#end

def standardize_probis_dict2(dict1, dict2):
    del_list = []
    for key in dict1:
        if dict1[key]== '-': #key is 3
            del_list.append(key)
    for key in del_list:
        del dict1[key]
    keysprobis = dict2.keys()
    keyscif = dict1.keys()
    final_dict = dict(zip(keysprobis, keyscif))
    return final_dict

def get_seq_index(pdb_id):
     protcif = PDBList().retrieve_pdb_file(pdb_id)
     pdb_info = MMCIF2Dict(protcif)
     protseq = pdb_info["_entity_poly.pdbx_seq_one_letter_code"]
     protseqfin = "".join(protseq).replace('\n', '')
     return protseqfin

full_list = []

#sort pdb, cen files
pdb_files = glob(CONFIG.probis_pdb_files)
pdb_files.sort()
cen_files = glob(CONFIG.probis_cen_files)
cen_files.sort()


for i,j in zip(pdb_files, cen_files):
     #create list for each binding site centroid, all stored in a big list
     cen_array = []
     with open(str(j), 'r') as file:
          contents = file.readlines()
          for line in contents:
               sample = str(line).split('\t')
               sample_arr = list(sample)
               del sample_arr[-1]
               cen_array.append(sample_arr)
          del cen_array[0]
     cen_array

#if __name__ == "__main__":
     #pdb_files = glob('/Users/akshome/Downloads/receptor*.pdb')
     #for fileName in pdb_files:
     structure_id = i.rsplit('/', 1)[1][:-4]
     structure = parser.get_structure(structure_id, i)[0]
     #print(list(structure))
     pdb_code = structure_id.rsplit('_', -1)[1]
     chain_id = structure_id.rsplit('_', -1)[2]
     #for chain in structure.get_chains():
     atoms = Bio.PDB.Selection.unfold_entities(structure, 'A')
     resislist = []
     probisres_dict = {}
     for chain in structure:
          for residue in chain:
               resi1 = residue.resname.title()
               resi11 = Bio.Data.IUPACData.protein_letters_3to1[resi1]
               resislist.append(resi11)
               probisres_dict[residue.get_full_id()[3][1]] = resi11
               resstr = "".join(resislist)
     ns = Bio.PDB.NeighborSearch(atoms)
     seq_list = []
     for item in cen_array:
          x,y,z = item[1:4]
          close_res = ns.search(np.array([x,y,z], float), float(item[4]), level='R')
          #res_list.extend(close_res)
          for residue in close_res:
               resseq = residue.get_full_id()[3][1]
               seq_list.append(resseq)
     seq_list = str(sorted(set(seq_list)))
     #res_list = sorted(set(res_list))
     prot_sequence = get_seq_index(pdb_code)
     seq_prot = Seq(prot_sequence)
     seq_res = Seq(resstr)
     alignments = aligner.align(seq_prot, seq_res)
     for alignment in alignments:
          probisalign = alignment[1]
     probis_dict = dict(zip(np.arange(1,len(probisalign)+1), probisalign))
     mapping_dict = standardize_probis_dict2(probis_dict, probisres_dict)
     seq_list2 =[ mapping_dict[i] for i in seq_list]
     seq_list2 = str(seq_list2)
     dflist = [pdb_code, chain_id, seq_list2]
     probis_tbl = add_data_to_example(con, probis_tbl, dflist)
     #full_list.append(dflist)
     #print(f"1: {pdb_code}  2: {chain_id}  3: {res_list} 4: {seq_list}")
     #print('End of file')
#probisdf = pd.DataFrame(
     #[i[0], i[1], i[3]] for i in full_list
#)
#probisdf = probisdf.rename({0: 'PDB ID', 1: 'Chain ID', 2: 'Sequence ID\'s of binding site residues'}, axis='columns')
print(pd.read_sql_query("SELECT * FROM example LIMIT 50", con))
     #print(type(pdb_code), type(chain.id), type(res_list))


