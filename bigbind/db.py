import os
import sqlite3

from config import CONFIG

def create_connection():
    """ create a database connection to the BigBind database """
    con = sqlite3.connect(CONFIG.db_file)
    # see https://stackoverflow.com/questions/9937713/does-sqlite3-not-support-foreign-key-constraints
    con.execute('PRAGMA foreign_keys = ON')
    return con

def create_tables():
    """ Creates an entirely new database, creating the tables
    from the create_tables.sql file. Returns a connection to the
    database. """

    with open("sql/create_tables.sql", "r") as f:
        sql = f.read()

    if os.path.exists(CONFIG.db_file):
        os.remove(CONFIG.db_file)
    con = create_connection()
    con.executescript(sql)

    return con