import csv
import MySQLdb

import rdkit
from rdkit import Chem
from rdkit.Chem import MACCSkeys

def get_cursor():
    conn = MySQLdb.connect(host = "localhost", user = "", passwd = "", db = "chembl_18")
    cursor = conn.cursor(cursorclass=MySQLdb.cursors.DictCursor)
    return cursor


def get_compounds():
	compounds = []
	cursor = get_cursor()
	with open("../source_data/chembl_gr.csv", "r") as inputfile:
		chembl_ids = [r["chembl_id"] for r in csv.DictReader(inputfile)]

	for c in chembl_ids:
		query = """SELECT chembl_id, canonical_smiles FROM compound_structures
					LEFT JOIN molecule_dictionary ON molecule_dictionary.molregno = compound_structures.molregno
					WHERE molecule_dictionary.chembl_id = '{}'""".format(c)
		cursor.execute(query)
		compounds.extend(cursor.fetchall())
	return compounds

def generate_maccs_keys(compounds):
	for i, c in enumerate(compounds):
		mol = Chem.MolFromSmiles(c["canonical_smiles"])
		compounds[i]["maccs"] = MACCSkeys.GenMACCSKeys(mol).ToBitString()
	return compounds

if __name__ == '__main__':
	compounds = get_compounds()
	compounds = generate_maccs_keys(compounds)
	print MACCSkeys.smartsPatts
	