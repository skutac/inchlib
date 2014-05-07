import csv

import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem


input_file = open("../source_data/chembl_era_ligands_descriptors.csv", "r")
reader = csv.DictReader(input_file, delimiter=",")

molid2smiles = {r["chembl_id"]: r["smiles"] for r in reader if r["standard_type"] == "Ki"}


for m in molid2smiles:
	mol = Chem.MolFromSmiles(molid2smiles[m])
	AllChem.Compute2DCoords(mol)
	Draw.MolToFile(mol,'../img/gr_molecules/{}.svg'.format(m))