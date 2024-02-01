import os
import sqlite3

DB_FILE = "data/bigbind.db"

def create_connection():
    """ create a database connection to the BigBind database """
    con = sqlite3.connect(DB_FILE)
    # see https://stackoverflow.com/questions/9937713/does-sqlite3-not-support-foreign-key-constraints
    con.execute('PRAGMA foreign_keys = ON')
    return con

def create_tables():
    """ Creates an entirely new database, creating the tables
    from the create_tables.sql file. Returns a connection to the
    database. """

    with open("create_tables.sql", "r") as f:
        sql = f.read()

    if os.path.exists(DB_FILE):
        os.remove(DB_FILE)
    con = create_connection()
    con.executescript(sql)

    return con