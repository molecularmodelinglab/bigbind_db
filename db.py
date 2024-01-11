import sqlite3

def create_connection():
    """ create a database connection to the BigBind database """
    return sqlite3.connect("data/bigbind.db")

def make_bigbind():
    con = create_connection()
    