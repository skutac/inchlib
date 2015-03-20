import csv, shutil

import rdkit
from rdkit import Chem
from rdkit.Chem import Draw, AllChem

def prepare_data_for_images_example():
	with open("static/source_data/steroids_fragments.csv", "r") as input_file:
		sid2fragments = {int(r[0]):r[1:] for r in csv.reader(input_file)}

	with open("static/source_data/steroids.csv", "r") as input_file:
		rows = [r for r in csv.reader(input_file)]

	filtered = [rows[0]]

	for r in rows[1:]:
		if int(r[0]) in sid2fragments:
			filtered.append(r)

	for r in filtered[1:]:
		for f in sid2fragments[int(r[0])]:
			if not "None" in f:
				f = f.strip("' ")
				shutil.copyfile("/home/ctibor/repositories/ChemGenDBD/chemgendbd/static/compounds/fragments/{}.svg".format(f), "static/img/fragments/{}.svg".format(f))

	with open("static/source_data/steroids_img_subset.csv", "w") as w:
		writer = csv.writer(w)
		writer.writerows(filtered)

def create_fragment_images():
	with open("static/source_data/steroids_fragments.csv", "r") as input_file:
		sid2fragments = {int(r[0]):r[1:] for r in csv.reader(input_file)}

	for sid, fs in sid2fragments.items():
		for f in fs:
			if f != "None":
				Draw.MolToFile(Chem.MolFromSmiles(f), "static/img/fragments/{}.svg".format(f), size=(100, 50))

def create_fragment_fps():
	with open("static/source_data/steroids_fragments.csv", "r") as input_file:
		sid2fragments = {int(r[0]):r[1:] for r in csv.reader(input_file)}

	sid2fps = {}

	for sid, fs in sid2fragments.items():
		sid2fps[sid] = []
		for f in fs:
			if f != "None":
				sid2fps[sid].append(get_morgan_fingerprint_for_smiles(f, radius=2, length=256))
			else:
				sid2fps[sid].append("".join(["0" for i in range(256)]))
	
	with open("static/source_data/steroids.csv", "r") as input_file:
		ids = [r[0] for r in csv.reader(input_file)]
	
	ids = set([int(v) for v in ids[1:]])

	with open("static/source_data/steroids_fragment_fps.csv", "w") as w:
		writer = csv.writer(w)
		writer.writerow(["id","Position 1","Position 2","Position 12","Position 13"])
		unique = set()

		for sid, fp in sid2fps.items():
			fp = "".join(fp)
			if not fp in unique and sid in ids:
				unique.add(fp)
				row = [sid]
				row.extend([i for i in fp])
				writer.writerow(row)

def get_morgan_fingerprint_for_smiles(smiles, radius=3, length=1024):
      mol = Chem.MolFromSmiles(smiles)
      fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=length)
      return fp.ToBitString()

if __name__ == '__main__':
	# prepare_data_for_images_example()
	# create_fragment_images()
	create_fragment_fps()