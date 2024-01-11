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

This query takes a while to execute! Make sure to save it to a csv in `temp_data` so you only need to execute this once.

This dataframe has many columns, but here are some important ones:
- `canonical_smiles`: The SMILES string representing the current molecule. Each _unique_ SMILES should correspond to a row in the `molecules` table. A first stab at creating the table will create a dataframe whose rows are the unique smiles, and with all the other rows (`has_conformer`, etc..) set to None.
- `protein_accession`: the UniProt ID of the protein target
- `standard_type`, `standard_relation`, `standard_units`, `standard_value`, and `pchembl_value`: the type of the activity measurement (IC50, Ki, etc), the relationship ('=', '>', or '<'; we only care about '='), the units (eg nM), the value, and the -log10 of the value (a more useful value in eg ML training)

Once we have the `molecules` table, the next steps are:

- Creating the `proteins` table. We want each row to correspond to a unqiue protein sequence. You can get the sequences from the uniprot IDs like [this](https://stackoverflow.com/questions/52569622/protein-sequence-from-uniprot-protein-id-python)
- Creating the `activities` table.
- Adding in the extra molecule information to `molecules`

# Loading ProBis data

Creating the actual `binding_sites` and associated files properly requires the `components` table, which we don't have yet. To prepare ourselves for when that occurs, we want a script that goes through all the pockets in ProBis and prints the pdb id, chain id, and residue ids of each site. Once we have all that info, creating the tables will be easy.

You can find all the ProBis data at `/proj/kpoplab/ProbisDock/` on longleaf

# Loading 
