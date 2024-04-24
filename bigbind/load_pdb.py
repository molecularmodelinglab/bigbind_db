import biotite.structure.io.pdbx as pdbx
import numpy as np
import pandas as pd
import os
import biotite
from rdkit import Chem
from rdkit.Chem import MolToSmiles, MolFromSmiles
from collections import Counter
from PyAstronomy import pyasl
an = pyasl.AtomicNo()
from collections import OrderedDict
import tqdm as tqdm
import sqlite3
from bigbind.db import create_connection
from bigbind.config import CONFIG
from chembl_structure_pipeline import standardizer
import duckdb
from bigbind.load_chembl import gotzinc
from rdkit import RDLogger

def chembl_smiles(bb_smiles):
    std_mol = standardizer.standardize_mol(MolFromSmiles(bb_smiles))
    return MolToSmiles(std_mol)

def proteins_from_cif(file_path):
    pdbx_file = pdbx.PDBxFile.read(file_path)
    poly = pdbx_file.get('entity_poly')

    prot_comps = []

    if type(poly['type']) == str:
        if poly['type'] != 'polydeoxyribonucleotide':
            for i in poly['pdbx_strand_id'].split(','):
                if 'X' in poly['pdbx_seq_one_letter_code_can']:
                    prot_comps.append([list(pdbx_file.keys())[0][0], i, -1, '', poly['pdbx_seq_one_letter_code'],'Protein'])
                else:
                    prot_comps.append([list(pdbx_file.keys())[0][0], i, -1,'', poly['pdbx_seq_one_letter_code_can'],'Protein'])
        else:
            for i in poly['pdbx_strand_id'].split(','):
                if 'X' in poly['pdbx_seq_one_letter_code_can']:
                    prot_comps.append([list(pdbx_file.keys())[0][0], i, -1, '', poly['pdbx_seq_one_letter_code'],'Nucleic Acid'])
                else:
                    prot_comps.append([list(pdbx_file.keys())[0][0], i, -1,'', poly['pdbx_seq_one_letter_code_can'],'Nucleic Acid'])
    else:
        for p in enumerate(poly['type']):
            if p[1] != 'polydeoxyribonucleotide':
                new_chains = poly['pdbx_strand_id'][p[0]].split(",")
                for n in new_chains:
                    if 'X' in poly['pdbx_seq_one_letter_code_can'][p[0]]:
                        prot_comps.append([list(pdbx_file.keys())[0][0], n, -1,'', poly['pdbx_seq_one_letter_code'][p[0]], 'Protein'])
                    else:
                        prot_comps.append([list(pdbx_file.keys())[0][0], n, -1,'', poly['pdbx_seq_one_letter_code_can'][p[0]], 'Protein'])
            else:
                new_chains = poly['pdbx_strand_id'][p[0]].split(",")
                for n in new_chains:
                    if 'X' in poly['pdbx_seq_one_letter_code_can'][p[0]]:
                        prot_comps.append([list(pdbx_file.keys())[0][0], n, -1,'', poly['pdbx_seq_one_letter_code'][p[0]], 'Nucleic Acid'])
                    else:
                        prot_comps.append([list(pdbx_file.keys())[0][0], n, -1,'', poly['pdbx_seq_one_letter_code_can'][p[0]], 'Nucleic Acid'])

    return prot_comps

def create_protcomps(dir_path):
    prots = []
    for x in os.listdir(dir_path):
        path1 = os.path.join(dir_path, x)
        singleprot = proteins_from_cif(path1)
        if singleprot != None:
            for m in singleprot:
                prots.append(m)
            
    df = pd.DataFrame(prots, columns = ['structure_id', 'chain', 'residue','smiles', 'sequence','type'])
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
                    mol_list.append([list(pdbx_file.keys())[0][0], chain, int(res), mol_from_stack(working_stack, working_bonds.as_array()), '','Ligand'])
        
        return mol_list

def create_ligcomps(dir_path):
    mols = []
    for x in os.listdir(dir_path):
        path1 = os.path.join(dir_path, x)
        singlelig = cif_to_mols(path1)
        if singlelig != None:
            for m in singlelig:
                mols.append(m)
            
    df = pd.DataFrame(mols, columns = ['structure_id', 'chain', 'residue', 'smiles','sequence','type'])
    return df


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
        structs_list.append([pdbid,pdbid, expmethod, resolution])
    structures = pd.DataFrame(structs_list, columns = ['id','pdb', 'type', 'resolution'])
    return structures

def create_tempcomps(dir_path):
    protcomps = create_protcomps(dir_path)
    ligcomps = create_ligcomps(dir_path)

    df = pd.concat([protcomps, ligcomps])
    df.insert(loc=0, column='id', value=np.arange(1, len(df)+1, dtype=int))
    
    return(df)

def create_components(dir_path, tempcomps):
    df = tempcomps.copy()
    df.drop('smiles', axis=1, inplace=True)
    df.drop('sequence', axis=1, inplace=True)
    return df

def create_ligand_components(tempcomps, con):

    df = duckdb.sql("SELECT id, smiles FROM tempcomps WHERE type = 'Ligand'").to_df()
    df.rename(columns={'id':'components_id'},inplace=True)
    df.insert(1, 'molecule_id', '')
    RDLogger.DisableLog('rdApp.*')
    for index, row in df.iterrows(): # room for potential efficiency gains w/ itertuples and such
        can_smile = chembl_smiles(row['smiles'])
        cur1 = con.cursor()
        cur1.execute("SELECT id FROM molecules WHERE smiles = ?",(can_smile,))
        result = cur1.fetchone()
        if result != None:
            df.at[index, 'molecule_id'] = result[0]
        else:
            sql = "INSERT INTO molecules(smiles,has_conformer,zinc_elements,num_components,molecular_weight,tanimoto_split,chembl_id) VALUES(?,?,?,?,?,?,?)"
            values = (can_smile, None, gotzinc(MolFromSmiles(can_smile)), 1, Chem.Descriptors.ExactMolWt(MolFromSmiles(can_smile)), '', '')
            cur2 = con.cursor()
            cur2.execute(sql, values)
            con.commit()
            cur3 = con.cursor()
            cur3.execute("SELECT id FROM molecules WHERE smiles =?",(can_smile,))
            result1 = cur3.fetchone()
            df.at[index, 'molecule_id'] = result1[0]
            
    return df

def create_protein_components(tempcomps, con):
    
    df = duckdb.sql("SELECT id, sequence FROM tempcomps WHERE type = 'Protein'").to_df()
    df.rename(columns={'id':'components_id'},inplace=True)
    df.insert(1, 'protein_id', '')
    df.insert(2, 'has_ptms', '')
    for index, row in df.iterrows():
        curp = con.cursor()
        curp.execute("SELECT id FROM proteins WHERE sequence = ?", (row['sequence'],))
        result = curp.fetchone()
        if result != None:
            df.at[index, 'protein_id'] = result[0]
        else:
            sql = "INSERT INTO proteins(sequence, uniprot, mutations, sequence_split) VALUES (?,?,?,?)"
            values = (row['sequence'], '','','')
            curp1 = con.cursor()
            curp1.execute(sql, values)
            con.commit()
            curp2 = con.cursor()
            curp2.execute("SELECT id FROM proteins WHERE sequence = ?", (row['sequence'],))
            result1 = curp2.fetchone()
            df.at[index, 'protein_id'] = result1[0]
    df.drop('sequence', axis=1, inplace=True)
    return df

def find_covalent_attachments(dir_path, con):
    #df = pd.DataFrame(columns=['component_1_id', 'component_2_id'])
    dnathings=['DA','DT','DG','DC','DU', 'TGP', '8OG']
    aminos = ['PTR','VAL', 'ILE', 'LEU', 'GLU', 'GLN', 'ASP', 'ASN', 'HIS', 'TRP', 'PHE', 'TYR', 'ARG', 'LYS', 'SER', 'THR', 'MET', 'ALA', 'GLY', 'PRO', 'CYS']
    for x in os.listdir(dir_path):
        filepath = os.path.join(dir_path, x)
        pdbx_file = pdbx.PDBxFile.read(filepath)
        pdbid = pdbx_file.get('exptl')['entry_id']
        if 'struct_conn' in pdbx_file.keys():
            d = pdbx_file.get('struct_conn')
            for y in enumerate(d['conn_type_id']):
                if y[1] == 'covale':
                    c1chain = d['ptnr1_auth_asym_id'][y[0]]
                    c1id = d['ptnr1_auth_seq_id'][y[0]]
                    c1res =d['ptnr1_auth_comp_id'][y[0]]
                    c2chain = d['ptnr2_auth_asym_id'][y[0]]
                    c2id = d['ptnr2_auth_seq_id'][y[0]]
                    c2res =d['ptnr2_auth_comp_id'][y[0]]
                    
                    curc1 = con.cursor()
                    if c1res in aminos:
                        c1type = "Protein"
                        curc1.execute('SELECT id FROM components WHERE structure_id = ? AND chain = ? AND type = ?', (pdbid, c1chain, c1type))
                        result1 = curc1.fetchone()[0]
                    elif c1res in dnathings:
                        c1type = "Nucleic Acid"
                        curc1.execute('SELECT id FROM components WHERE structure_id = ? AND chain = ? AND type = ?', (pdbid, c1chain, c1type))
                        result1 = curc1.fetchone()[0]
                    else:
                        c1type = "Ligand"
                        curc1.execute('SELECT id FROM components WHERE structure_id = ? AND chain = ? AND residue = ? AND type = ?', (pdbid, c1chain, c1id, c1type))
                        print(pdbid, c1chain, c1res, c1id, c1type)
                        result1 = curc1.fetchone()[0]
                    
                    curc2 = con.cursor()
                    if c2res in aminos:
                        c2type = "Protein"
                        curc2.execute('SELECT id FROM components WHERE structure_id = ? AND chain = ? AND type = ?', (pdbid, c2chain, c2type))
                        result2 = curc2.fetchone()[0]
                    elif c2res in dnathings:
                        c2type = "Nucleic Acid"
                        curc2.execute('SELECT id FROM components WHERE structure_id = ? AND chain = ? AND type = ?', (pdbid, c2chain, c2type))
                        result2 = curc2.fetchone()[0]
                    else:
                        c2type = "Ligand"
                        curc2.execute('SELECT id FROM components WHERE structure_id = ? AND chain = ? AND residue = ? AND type = ?', (pdbid, c2chain, c2id, c2type))
                        print(pdbid, c2chain, c2res, c2id, c2type)
                        result2 = curc2.fetchone()[0]
                        
                    
                    curs = con.cursor()
                    sql = "INSERT INTO covalent_attachments(component_1_id, component_2_id) VALUES (?, ?)"
                    values = (int(result1), int(result2))
                    print(values)
                    curs.execute(sql, values)
                    con.commit()       


def load_pdb():

    pdbs = "C:\\Users\\anees\\Downloads\\randomfolder1" # change to whatever path of cifs is later

    con = create_connection()
    cur = con.cursor()

    

    structures = create_structures(pdbs)
    temp_comps = create_tempcomps(pdbs)
    components = create_components(pdbs, temp_comps)

    ligand_components = create_ligand_components(temp_comps, con)
    protein_components = create_protein_components(temp_comps, con)

    

    structures.to_sql(con=con, name='structures', schema='SCHEMA', index=False, if_exists='replace')
    components.to_sql(con=con, name='components', schema='SCHEMA', index=False, if_exists='replace')
    protein_components.to_sql(con=con, name='protein_components', schema='SCHEMA', index=False, if_exists='replace')
    ligand_components.to_sql(con=con, name='ligand_components', schema='SCHEMA', index=False, if_exists='replace')

    #find_covalent_attachments(pdbs, con)

    #cur.execute("SELECT * FROM components")
    #result = cur.fetchall()
    #for row in result: 
    #    print(row) 
    #    print("\n") 


if __name__ == "__main__":
    load_pdb()
    