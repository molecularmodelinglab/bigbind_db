import biotite.structure.io.pdbx as pdbx
import numpy as np
import pandas as pd
import os
import biotite
from rdkit import Chem
from rdkit.Chem import MolToSmiles
from collections import Counter
from PyAstronomy import pyasl
an = pyasl.AtomicNo()
from collections import OrderedDict
import tqdm as tqdm
import sqlite3
from bigbind.db import create_connection
from bigbind.config import CONFIG




def proteins_from_cif(file_path):
    pdbx_file = pdbx.PDBxFile.read(file_path)
    poly = pdbx_file.get('entity_poly')

    prot_comps = []

    if type(poly['type']) == str:
        if poly['type'] != 'polydeoxyribonucleotide':
            for i in poly['pdbx_strand_id'].split(','):
                if 'X' in poly['pdbx_seq_one_letter_code_can']:
                    prot_comps.append([list(pdbx_file.keys())[0][0], poly['pdbx_seq_one_letter_code'], i])
                else:
                    prot_comps.append([list(pdbx_file.keys())[0][0], poly['pdbx_seq_one_letter_code_can'], i])
    else:
        for p in enumerate(poly['type']):
            if p[1] != 'polydeoxyribonucleotide':
                new_chains = poly['pdbx_strand_id'][p[0]].split(",")
                for n in new_chains:
                    if 'X' in poly['pdbx_seq_one_letter_code_can'][p[0]]:
                        prot_comps.append([list(pdbx_file.keys())[0][0], poly['pdbx_seq_one_letter_code'][p[0]], n])
                    else:
                        prot_comps.append([list(pdbx_file.keys())[0][0], poly['pdbx_seq_one_letter_code_can'][p[0]], n])

    return prot_comps

def create_protcomps(dir_path):
    prots = []
    for x in os.listdir(dir_path):
        path1 = os.path.join(dir_path, x)
        singleprot = proteins_from_cif(path1)
        if singleprot != None:
            for m in singleprot:
                prots.append(m)
            
    df = pd.DataFrame(prots, columns = ['Entry', 'Sequence', 'Chain'])
    return df

def add_formal_charges(m):
    # lol
    m.UpdatePropertyCache(strict=False)
    for at in m.GetAtoms():
        if at.GetAtomicNum() == 7 and at.GetExplicitValence()==4 and at.GetFormalCharge()==0:
            at.SetFormalCharge(1)
        if at.GetAtomicNum() == 8 and at.GetExplicitValence()==3:
            at.SetFormalCharge(1)
        if at.GetAtomicNum() == 8 and at.GetExplicitValence()==1 and at.GetFormalCharge()==0:
            at.SetFormalCharge(-1)
        if at.GetAtomicNum() == 5 and at.GetExplicitValence()==4 and at.GetFormalCharge()==0:
            at.SetFormalCharge(-1)
        if at.GetAtomicNum() == 8 and at.GetExplicitValence()==2 and at.GetFormalCharge()!=0:
            at.SetFormalCharge(0)


# takes stack and makes rdkit mol object
def mol_from_stack (stack, bonds):
    mol = Chem.RWMol()
    name2idx = {}
    b = 0
    for atom in stack:
        atom1 = Chem.Atom(an.getAtomicNo(str(atom.element).title()))
        atom1.SetFormalCharge(int(atom.charge))
        i = mol.AddAtom(atom1)
        name2idx[b] = i
        b += 1
    for bond in bonds:
        i1 = name2idx[bond[0]]
        i2 = name2idx[bond[1]]
        order = {
            1: Chem.BondType.SINGLE,
            2: Chem.BondType.DOUBLE,
            3: Chem.BondType.TRIPLE,
            5: Chem.BondType.SINGLE,
            6: Chem.BondType.DOUBLE}[bond[2]]
        mol.AddBond(i1, i2, order)
    # use openbabel to protonate everything at pH seven, THEN add formal charges
    add_formal_charges(mol)
    Chem.SanitizeMol(mol)
    smile = MolToSmiles(mol)
    return smile



def cif_to_mols(file_path):
    # load file
    pdbx_file = pdbx.PDBxFile.read(file_path)
    stack = pdbx.get_structure(pdbx_file, extra_fields=["charge"], model=1)
    bonds = biotite.structure.connect_via_residue_names(stack)

    # nonpolymer chain info for parsing
    if 'pdbx_nonpoly_scheme' in pdbx_file.keys():
        nonpoly = pdbx_file.get('pdbx_nonpoly_scheme')
    
       
        cres = Counter(nonpoly['mon_id'])

        # i think this is inefficient and i can just make a np array w zeros atp
        boolList = []
        important_list = []
        for atom in stack:
            if any(item == atom.res_name for item in list(cres.keys())) and atom.res_name != 'HOH':
                boolList.append(1)
                important_list.append([atom.chain_id, atom.res_id, atom.res_name])
            else:
                boolList.append(0)
        # dunno why i did this
        boolArray = np.array(boolList, dtype=bool)

        # stack of only small molecules
        new_stack = stack[boolArray]

        
        # this is the current return, but will change when I figure out how exactly we want to populate the df
        mol_list = []

        for chain in list(OrderedDict.fromkeys(l[0] for l in important_list)):
            for res in list(OrderedDict.fromkeys(x[1] for x in important_list)):
                new_slice = []
                for atom in new_stack:
                    if atom.res_id == res and atom.chain_id == chain:
                        new_slice.append(1)
                    else:
                        new_slice.append(0)
                slice_array = np.array(new_slice, dtype=bool)
                working_stack = new_stack[slice_array]
                if working_stack.shape != (0,):
                    resname = working_stack[0].res_name
                    working_bonds = biotite.structure.connect_via_residue_names(working_stack)
                    mol_list.append([list(pdbx_file.keys())[0][0], chain, resname, res, mol_from_stack(working_stack, working_bonds.as_array())])
        
        return mol_list

def create_ligcomps(dir_path):
    mols = []
    for x in os.listdir(dir_path):
        path1 = os.path.join(dir_path, x)
        singlelig = cif_to_mols(path1)
        if singlelig != None:
            for m in singlelig:
                mols.append(m)
            
    df = pd.DataFrame(mols, columns = ['Entry', 'Chain', 'Residue','Res_ID', 'smiles'])
    return df

#def create_temp_comps(dir_path):

def create_structures(dir_path):
    structs_list = []
    for x in os.listdir(dir_path):
        path1 = os.path.join(dir_path, x)
        pdbx_file = pdbx.PDBxFile.read(path1)
        pdbid = pdbx_file.get('exptl')['entry_id']
        expmethod = pdbx_file.get('exptl')['method']
        try:
            resolution = float(pdbx_file.get('refine')['ls_d_res_high'])
        except:
            resolution = float(pdbx_file.get('em_3d_reconstruction')['resolution'])
        structs_list.append([pdbid, expmethod, resolution])
    structures = pd.DataFrame(structs_list, columns = ['pdb', 'type', 'resolution'])
    return structures



def load_pdb():

    pdbs = "C:\\Users\\anees\\Downloads\\randomfolder1" # change to whatever path of cifs is later

    #protein_components = create_protcomps(pdbs)

    #ligand_components = create_ligcomps(pdbs)

    structures = create_structures(pdbs)

    con = create_connection()
    
    structures.to_sql(con=con, name='structures', schema='SCHEMA', index=True, index_label='id', if_exists='replace')
    #protein_components.to_sql(con=con, name='protein_components', schema='SCHEMA', index=False, if_exists='append')
    #ligand_components.to_sql(con=con, name='ligand_components', schema='SCHEMA', index=False, if_exists='append')

if __name__ == "__main__":
    load_pdb()
    