import csv
import MySQLdb

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Scaffolds import MurckoScaffold

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

def get_scaffolds(compounds):
	for i, c in enumerate(compounds):
		mol = Chem.MolFromSmiles(c["canonical_smiles"])
		core = MurckoScaffold.GetScaffoldForMol(mol)
		compounds[i]["scaffold"] = Chem.MolToSmiles(core)
		compounds[i]["generic_scaffold"] = Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(core))
	return compounds

def create_scaffold_img(scaffold2num):
	for m in scaffold2num:
		mol = Chem.MolFromSmiles(m)
		AllChem.Compute2DCoords(mol)
		Draw.MolToFile(mol,'../img/gr_scaffolds/{}.svg'.format(scaffold2num[m]))

def export_compound2scaffold(compounds, scaffold2num):
	with open("../source_data/compound2scaffold.csv", "w") as output:
		output.write("chembl_id,scaffold\n")
		for c in compounds:
			output.write(",".join([c["chembl_id"], str(scaffold2num[c["generic_scaffold"]])]))
			output.write("\n")
	return


if __name__ == '__main__':
	compounds = get_compounds()
	scaffolds = get_scaffolds(compounds)
	unique_scaffolds = list(set([c["generic_scaffold"] for c in compounds]))
	unique_scaffolds.sort()
	scaffold2num = {s:i for i,s in enumerate(unique_scaffolds)}

	create_scaffold_img(scaffold2num)
	export_compound2scaffold(compounds, scaffold2num)