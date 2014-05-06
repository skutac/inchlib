import csv

import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem


input_file = open("../source_data/chembl_GR_ligands_descriptors2.csv", "r")
reader = csv.DictReader(input_file, delimiter=",")

molid2smiles = {r["chembl_id"]: r["smiles"] for r in reader if r["standard_type"] == "IC50"}


for m in molid2smiles:
	mol = Chem.MolFromSmiles(molid2smiles[m])
	AllChem.Compute2DCoords(mol)
	Draw.MolToFile(mol,'../img/gr_molecules/{}.svg'.format(m))