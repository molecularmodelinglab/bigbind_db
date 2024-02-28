import biotite.structure.io.pdbx as pdbx
import numpy as np
import pandas as pd
import os

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

if __name__ == '__main__':
    print('hii')
    dir_path = "C:\\Users\\anees\\Downloads\\randomfolder1" # insert file path here
    prots = []
    for x in os.listdir(dir_path):
        path1 = os.path.join(dir_path, x)
        singleprot = proteins_from_cif(path1)
        if singleprot != None:
            for m in singleprot:
                prots.append(m)
            
    df = pd.DataFrame(prots, columns = ['Entry', 'Sequence', 'Chain'])
    df.to_csv('protein_comp.csv')
    print(df)