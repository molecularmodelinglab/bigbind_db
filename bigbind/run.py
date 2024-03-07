from prefect import flow

from bigbind.load_chembl import load_chembl
from bigbind.db import create_tables

@flow
def make_bigbind():
    create_tables()
    chembl_future = load_chembl()

if __name__ == "__main__":
    make_bigbind()