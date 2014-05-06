import csv
import MySQLdb
import numpy

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

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

def filter_gr_data(filters):
	with open("../source_data/chembl_GR_ligands_descriptors2.csv", "r") as gr_data:
		rows = [r for r in csv.DictReader(gr_data)]

	ids = [r["chembl_id"] for r in rows]

	for f in filters:
		rows = [r for r in rows if r[f] == filters[f]]


	mol2values = {r["chembl_id"]:[] for r in rows}

	for r in rows:
		mol2values[r["chembl_id"]].append(float(r["nm"]))

	rows = {r["chembl_id"]:r for r in rows}
	final = []

	for m in mol2values:
		dev = numpy.std(mol2values[m])
		if dev < 100:
			rows[m]["nm"] = numpy.median(mol2values[m])
			final.append(rows[m])

	return final

def export_gr_data(data, id_field, metadata_field):
	with open("../source_data/chembl_gr_2.csv", "w") as output:
		header = ["id"]
		header.extend(FIELDS)
		output.write(",".join(header))
		output.write("\n")

		for r in data:
			row = [str(r[id_field])]

			for f in FIELDS:
				row.append(str(round(float(r[f]), 2)))

			output.write(",".join(row))
			output.write("\n")

	with open("../source_data/chembl_gr_metadata_2.csv", "w") as output:
		header = ["id", metadata_field]
		output.write(",".join(header))
		output.write("\n")

		for r in data:
			output.write(",".join([str(r[id_field]), str(r[metadata_field])]))
			output.write("\n")


if __name__ == '__main__':
	FIELDS = ["SlogP", "NumHBD", "NumHBA", "ExactMW", "NumRotatableBonds"]
	FIELDS = ["SlogP", "SMR", "LabuteASA", "TPSA", "ExactMW", "NumRotatableBonds", "NumHBD", "NumHBA", "NumAmideBonds", "NumHeteroAtoms", "NumHeavyAtoms", "NumAtoms", "NumRings", "NumAromaticRings", "NumSaturatedRings", "NumAliphaticRings", "NumAromaticHeterocycles", "NumSaturatedHeterocycles", "NumAliphaticHeterocycles", "NumAromaticCarbocycles", "NumSaturatedCarbocycles", "NumAliphaticCarbocycles"]
	# FIELDS = ["MQN{}".format(i) for i in range(1, 42)]
	data = filter_gr_data({"standard_type": "IC50", "standard_relation": "=", "confidence_score": "9"})
	# data = filter_gr_data({"standard_type": "Potency", "confidence_score": "9"})
	export_gr_data(data, "chembl_id", "nm")
	# compounds = get_compounds()
	# scaffolds = get_scaffolds(compounds)
	# unique_scaffolds = list(set([c["generic_scaffold"] for c in compounds]))
	# unique_scaffolds.sort()
	# scaffold2num = {s:i for i,s in enumerate(unique_scaffolds)}

	# create_scaffold_img(scaffold2num)
	# export_compound2scaffold(compounds, scaffold2num)