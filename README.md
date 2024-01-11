# The BigBind Database

# Schema
We're gonna follow [this schema](https://dbdiagram.io/d/BigBind-65845d7989dea627995d0d84), modifiying it as needed. The idea is that making each of these tables is a small, contained step.

# Loading in data from ChEMBL

You can download the full chembl database as an SQLite file [here](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_33_schema.png). The schema is [here](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_33_schema.png). It's complicated!

Following the old bigbind code [here](https://github.com/molecularmodelinglab/bigbind/blob/7c67d75d7f4b0d0f74df5faa62372afc56d34c3d/bigbind/bigbind.py#L214), you can get a big pandas dataframe of all the activities (with some filtering) with the following SQL query:
```sql

SELECT md.chembl_id AS compound_chembl_id,
    cs.canonical_smiles,
    act.standard_type,
    act.standard_relation,
    act.standard_value,
    act.standard_units,
    act.pchembl_value,
    act.potential_duplicate,
    COALESCE(act.data_validity_comment, 'valid') as data_validity_comment,
    a.confidence_score,
    td.chembl_id AS target_chembl_id,
    td.target_type,
    c.accession as protein_accession,
    a.chembl_id as assay_id
    FROM target_dictionary td
    JOIN assays a ON td.tid = a.tid
    JOIN activities act ON a.assay_id = act.assay_id
    JOIN molecule_dictionary md ON act.molregno = md.molregno
    JOIN compound_structures cs ON md.molregno = cs.molregno
    JOIN target_type tt ON td.target_type = tt.target_type
    JOIN target_components tc ON td.tid = tc.tid
    JOIN component_sequences c ON tc.component_id = c.component_id
    AND tt.parent_type  = 'PROTEIN' 
    AND act.pchembl_value IS NOT NULL

```