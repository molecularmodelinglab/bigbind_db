# Example on how loading data into BigBind should work
# Basically put functions into prefect tasks, especially if
# they can be called concurrently

from prefect import task, flow
from bigbind.config import CONFIG

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