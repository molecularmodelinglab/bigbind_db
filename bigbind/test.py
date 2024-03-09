import numpy as np
import pandas as pd
import warnings
from bigbind.probis_tables import add_example_table, add_data_to_example
from bigbind.db import create_connection

con = create_connection()
test_tbl = add_example_table(con)
testlist = [1, 2, str([1, 2, 3])]
testlist2 = [1, 8, str([3, 4, 5])]
test_tbl = add_data_to_example(con, test_tbl, testlist)
test_tbl = add_data_to_example(con, test_tbl, testlist2)

print(pd.read_sql_query("SELECT * FROM example", con))
print(pd.read_sql_query("SELECT group_concat(residues) FROM example WHERE pdbID := pdb GROUP BY pdbID, {'pdb' : 101m}"), con)
