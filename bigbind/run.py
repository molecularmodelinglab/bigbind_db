from prefect import flow

from bigbind.load_chembl import load_chembl
from bigbind.db import create_tables
from bigbind.tanimoto_matrix import load_tanimoto_matrix

# @flow
def make_bigbind():
    create_tables()
    load_chembl()
    load_tanimoto_matrix()

if __name__ == "__main__":
    make_bigbind()