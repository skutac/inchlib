import csv

from propy import PyPro
from propy.GetProteinFromUniprot import GetProteinSequence

char2aa = {}

def get_sequence(pdbid):
	return GetProteinSequence(pdbid)

def calculate_descriptors(sequence):
	DesObject=PyPro.GetProDes(sequence) ##construct a GetProDes object
	descriptors = DesObject.GetAAComp() ##calculate 20 amino acid composition descriptors


def read_pdb_report(filename):
	with open(filename, "r") as input_file:
		pdbs = [r for r in csv.DictReader(input_file, delimiter=",")]

	return pdbs

if __name__ == '__main__':
	pdbs = read_pdb_report("../source_data/proteins_report.csv")
	
	for p in pdbs[:1]:
		print get_sequence(p["PDB ID"])
		# print calculate_descriptors(p["Sequence"])
		# p["Sequence"] = get_sequence(p["PDB ID"])

	# print pdbs