CREATE TABLE `proteins` (
  `id` integer PRIMARY KEY,
  `sequence` varchar,
  `uniprot` varchar,
  `mutations` varchar,
  `sequence_split` varchar,
  CONSTRAINT `sequence_unique` UNIQUE (`sequence`)
);

CREATE TABLE `molecules` (
  `id` integer PRIMARY KEY,
  `smiles` varchar,
  `has_conformer` boolean,
  `zinc_elements` boolean,
  `num_components` integer,
  `molecular_weight` integer,
  `tanimoto_split` varchar,
  `chembl_id` varchar,
  CONSTRAINT `smiles_unique` UNIQUE (`smiles`)
);

CREATE TABLE `activities` (
  `protein_id` integer,
  `ligand_id` integer,
  `type` varchar,
  `activity` float,
  `standard_relation` varchar,
  `standard_type` varchar,
  `standard_value` float, 
  `standard_units` varchar
  -- FOREIGN KEY (`ligand_id`) REFERENCES `molecules` (`id`),
  -- FOREIGN KEY (`protein_id`) REFERENCES `proteins` (`id`)
);

CREATE TABLE `structures` (
  `id` varchar PRIMARY KEY,
  `pdb` varchar,
  `type` varchar,
  `resolution` float,
  CONSTRAINT `pdb_unique` UNIQUE (`pdb`)
);

CREATE TABLE `components` (
  `id` integer PRIMARY KEY,
  `structure_id` varchar,
  `chain` varchar,
  `residue` integer,
  `type` varchar,
  FOREIGN KEY (`structure_id`) REFERENCES `structures` (`id`),
  CONSTRAINT `chain_residue_unique` UNIQUE (`structure_id`, `chain`, `residue`)
);

CREATE TABLE `protein_components` (
  `component_id` integer,
  `protein_id` integer,
  `has_ptms` boolean,
  FOREIGN KEY (`component_id`) REFERENCES `components` (`id`),
  FOREIGN KEY (`protein_id`) REFERENCES `proteins` (`id`)
);

CREATE TABLE `ligand_components` (
  `component_id` integer,
  `molecule_id` integer,
  `smiles` varchar,
  FOREIGN KEY (`component_id`) REFERENCES `components` (`id`),
  FOREIGN KEY (`molecule_id`) REFERENCES `molecules` (`id`)
);

CREATE TABLE `covalent_attachments` (
  `component_1_id` integer,
  `component_2_id` integer,
  FOREIGN KEY (`component_1_id`) REFERENCES `components` (`id`),
  FOREIGN KEY (`component_2_id`) REFERENCES `components` (`id`)
);

CREATE TABLE `binding_sites` (
  `id` integer PRIMARY KEY,
  `consensus_site` integer,
  `pocket_tm_split` varchar,
  FOREIGN KEY (`consensus_site`) REFERENCES `consensus_sites` (`id`)
);

CREATE TABLE `interaction_sites` (
  `id` integer PRIMARY KEY,
  `binding_site_id` integer,
  `includes_crystal_mate` boolean,
  FOREIGN KEY (`binding_site_id`) REFERENCES `binding_sites` (`id`)
);

CREATE TABLE `consensus_sites` (
  `id` integer PRIMARY KEY
);

CREATE TABLE `binding_site_components` (
  `binding_site_id` integer,
  `component_id` integer,
  `residues` varchar,
  FOREIGN KEY (`binding_site_id`) REFERENCES `binding_sites` (`id`),
  FOREIGN KEY (`component_id`) REFERENCES `components` (`id`)
);

CREATE TABLE `interaction_site_components` (
  `interaction_site_id` integer,
  `component_id` integer,
  `residues` varchar,
  FOREIGN KEY (`interaction_site_id`) REFERENCES `interaction_sites` (`id`),
  FOREIGN KEY (`component_id`) REFERENCES `components` (`id`)
);

CREATE TABLE `molecule_similarity` (
  `mol_1_id` integer,
  `mol_2_id` integer,
  `similarity` float,
  FOREIGN KEY (`mol_1_id`) REFERENCES `molecules` (`id`),
  FOREIGN KEY (`mol_2_id`) REFERENCES `molecules` (`id`)
);

CREATE TABLE `protein_similarity` (
  `prot_1_id` integer,
  `prot_2_id` integer,
  `similarity` float,
  FOREIGN KEY (`prot_1_id`) REFERENCES `proteins` (`id`),
  FOREIGN KEY (`prot_2_id`) REFERENCES `proteins` (`id`)
);

CREATE TABLE `binding_site_similarity` (
  `site_1_id` integer,
  `site_2_id` integer,
  `similarity` float,
  FOREIGN KEY (`site_1_id`) REFERENCES `binding_sites` (`id`),
  FOREIGN KEY (`site_2_id`) REFERENCES `binding_sites` (`id`)
);
