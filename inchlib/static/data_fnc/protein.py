import csv

from propy import PyPro
from propy.GetProteinFromUniprot import GetProteinSequence

char2aa = {
	"A":"Ala",
	"R":"Arg",
	"N":"Asn",
	"D":"Asp",
	"C":"Cys",
	"E":"Glu",
	"Q":"Gln",
	"G":"Gly",
	"H":"His",
	"I":"Ile",
	"L":"Leu",
	"K":"Lys",
	"M":"Met",
	"F":"Phe",
	"P":"Pro",
	"S":"Ser",
	"T":"Thr",
	"W":"Trp",
	"Y":"Tyr",
	"V":"Val"
}

def get_sequence(pdbid):
	return GetProteinSequence("Q99527")

def calculate_descriptors(sequence):
	DesObject=PyPro.GetProDes(sequence) ##construct a GetProDes object
	descriptors = DesObject.GetAAComp() ##calculate 20 amino acid composition descriptors
	return descriptors


def read_pdb_report(filename):
	with open(filename, "r") as input_file:
		pdbs = [r for r in csv.DictReader(input_file, delimiter=",")]

	return pdbs

def translate_descs(descs):
	return {char2aa[d]:descs[d] for d in descs}

if __name__ == '__main__':
	pdbs = read_pdb_report("../source_data/proteins_report.csv")
	pdb2desc = {}
	
	for p in pdbs:
		seq = []
		rows = p["Sequence"].split("\n")
		for r in rows:
			row_split = r.split(" ")
			for i in row_split:
				if not "1" in i and "" != i:
					seq.append(i.split("\n")[0])
		seq = "".join(seq)
		descs = calculate_descriptors(seq) 
		pdb2desc[p["PDB ID"]] = translate_descs(descs)


	keys = char2aa.values()
	header = ["PDB_ID"]
	header.extend(keys)
	
	with open("../source_data/proteins.csv", "w") as output:
		output.write(",".join(header))
		output.write("\n")

		for p in pdb2desc:
			row = []
			row.append(p)
			for k in keys:
				row.append(str(pdb2desc[p][k]))

			output.write(",".join(row))
			output.write("\n")