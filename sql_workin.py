# imports :)
import subprocess
import os
from glob import glob
import sqlite3
import pandas as pd
from tqdm import tqdm
import urllib.request
from pyfaidx import Fasta
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.ipython_3d = True
import py3Dmol
from rdkit.Chem import rdDepictor
from rdkit.Chem import rdDistGeom
import rdkit
from rdkit.Chem import AllChem
from pathlib import Path
from prefect import flow, task, get_run_logger
from bigbind.db import create_connection
from rdkit.Chem import Descriptors



if __name__ == "__main__":
    df = pd.read_csv("molecules.csv")
    con = create_connection()
    df.to_sql(con=con, name='molecules', schema='SCHEMA', index=False, if_exists='append')