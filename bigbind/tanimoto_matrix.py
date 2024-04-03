from array import array
from collections import defaultdict
from functools import partial
from multiprocessing import Pool, shared_memory, cpu_count
import random
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import pandas as pd
from scipy import sparse
import scipy
from tqdm import tqdm
import yaml
from bigbind.db import create_connection
from prefect import task, flow

def get_morgan_fps(cfg, df, radius=3, bits=2048):
    fps = []
    for smi in tqdm(df.smiles):
        mol = Chem.MolFromSmiles(smi)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, useChirality=True, radius=radius, nBits=bits)
        fps.append(np.array(fp, dtype=bool))
    fps = np.asarray(fps)
    return fps

def get_fp(smi, radius=4, bits=2048):
    mol = Chem.MolFromSmiles(smi)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, useChirality=True, radius=radius, nBits=bits)
    return np.array(fp, dtype=bool)

# like get_morgan_fps, but parallelized
MORGAN_FP_CPUS = cpu_count()//2
@task(timeout_seconds=1)
def get_morgan_fps_parallel(cfg, df):
    smi_list = list(df.smiles.unique())
    fps = []
    with Pool(MORGAN_FP_CPUS) as p:
        fps = list(tqdm(p.imap(get_fp, smi_list), total=len(smi_list)))
    fps = np.asarray(fps)
    return smi_list, fps

@task(timeout_seconds=1)
def get_morgan_fps_parallel_struct(cfg, df):
    smi_list = list(df.lig_smiles.unique())
    fps = []
    with Pool(MORGAN_FP_CPUS) as p:
        fps = list(tqdm(p.imap(get_fp, smi_list), total=len(smi_list)))
    fps = np.asarray(fps)
    return smi_list, fps

def batch_tanimoto(fp, fps):
    inter = np.logical_and(fp, fps)
    union = np.logical_or(fp, fps)
    sims = inter.sum(-1)/union.sum(-1)
    return sims

def get_bytes(a):
    return a.data.nbytes + a.row.nbytes + a.col.nbytes

def batch_tanimoto_faster(fp_shape, fp_shm_name, fp_sum_shape, fp_sum_shm_name, idx):
    fp_shm = shared_memory.SharedMemory(name=fp_shm_name)
    fps = np.ndarray(fp_shape, dtype=bool, buffer=fp_shm.buf)
    fp_sum_shm = shared_memory.SharedMemory(name=fp_sum_shm_name)
    fp_sum = np.ndarray(fp_sum_shape, dtype=int, buffer=fp_sum_shm.buf)

    fp = fps[idx]
    inter = np.logical_and(fp, fps)
    inter_sum = inter.sum(-1)
    sims = inter_sum/(fp_sum + fp.sum() - inter_sum)

    # print(idx, (sims < 0.2).sum()/len(sims))
    
    sims[sims < 0.4] = 0.0
    ssim = sparse.coo_matrix(sims)

    # print(ssim.data.shape, get_bytes(ssim))

    return ssim

TANIMOTO_CPUS = cpu_count()//2
def get_tanimoto_matrix_impl(fps):
    try:
        fp_sum = fps.sum(-1)

        fp_shm = shared_memory.SharedMemory(create=True, size=fps.nbytes)
        fps_shared = np.ndarray(fps.shape, dtype=bool, buffer=fp_shm.buf)
        fps_shared[:] = fps[:]
        
        fp_sum_shm = shared_memory.SharedMemory(create=True, size=fp_sum.nbytes)
        fp_sum_shared = np.ndarray(fp_sum.shape, dtype=int, buffer=fp_sum_shm.buf)
        fp_sum_shared[:] = fp_sum[:]

        sim_func = partial(batch_tanimoto_faster, fps.shape, fp_shm.name, fp_sum.shape, fp_sum_shm.name)
        with Pool(TANIMOTO_CPUS) as p:
            cols = list(tqdm(p.imap(sim_func, range(len(fps))), total=len(fps)))
        return sparse.vstack(cols)
    finally:
        fp_shm.close()
        fp_shm.unlink()
        fp_sum_shm.close()
        fp_sum_shm.unlink()

@task(timeout_seconds=24)
def get_tanimoto_matrix(cfg, fps):
    return get_tanimoto_matrix_impl(fps)

@task(timeout_seconds=24)
def get_tanimoto_matrix_struct(cfg, fps):
    return get_tanimoto_matrix_impl(fps)


@task()
def get_full_tanimoto_matrix(cfg, activities, smi_list, lig_sim_mat):
    """ Rejiggers the tanimoto matrix so that it's indexed by
    the indices of the activities dataframe, not the unique lig
    smiles index """
        
    smi2idx = { smi: idx for idx, smi in enumerate(smi_list) }
    idx2act_idx= defaultdict(set)
    for act_idx, smi in enumerate(activities.lig_smiles):
        idx = smi2idx[smi]
        idx2act_idx[idx].add(act_idx)

    new_row = array('I')
    new_col = array('I')
    new_data = array('f')
    for i, j, data in zip(tqdm(lig_sim_mat.row), lig_sim_mat.col, lig_sim_mat.data):
        for i2 in idx2act_idx[i]:
            for j2 in idx2act_idx[j]:
                new_row.append(i)
                new_col.append(j)
                new_data.append(data)

    new_row = np.array(new_row)
    new_col = np.array(new_col)
    new_data = np.array(new_data)

    return scipy.sparse.coo_array((new_data, (new_row, new_col)), shape=lig_sim_mat.shape, copy=False)

#gets matrix and connects to the sql db and outputs it there
def matrix_to_sql(matrix):
    con = create_connection
    print("populating df")
    nonzeroinx = matrix.nonzero()
    mol1id = nonzeroinx[0]
    mol2id = nonzeroinx[1]
    arr = matrix.toarray()
    val = []
    for i in tqdm(range(len(mol1id))):
        val.append(arr[mol1id[i]][mol2id[i]])
    
    
    
    df = pd.DataFrame()
    df["mol_1_id"] = mol1id
    df["mol_2_id"] = mol2id
    df["similarity"] = val
    df.to_sql(con=con, name='molecule_similarity', schema='SCHEMA', index=False, if_exists='append')


#main
@flow
def load_tanimoto_matrix():
    #set seed
    SEED = 49
    random.seed(SEED)
    np.random.seed(SEED)
    
    
    #get configs
    with open("configs/default.yml", "r") as f:
        cfg = yaml.safe_load(f)
    
    #tmp csv 
    #TODO CHANGE THIS TO ACTUALLY GETTING A QUERY FOR THE TABLE ON SQL
    df = pd.read_csv(f"{cfg['data_dir']}/csv_of_df/csv_files/csv_files/molecules.csv")
    #actual jazz
    fps = get_morgan_fps_parallel(cfg, df)
    tm = get_tanimoto_matrix(cfg, fps[1])
    
    matrix_to_sql(tm)

# SEED = 49
# random.seed(SEED)
# np.random.seed(SEED)
if __name__ == "__main__":
    # with open("configs/default.yml", "r") as f:
    #     cfg = yaml.safe_load(f)
    load_tanimoto_matrix()