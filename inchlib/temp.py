import csv, shutil

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

if __name__ == '__main__':
	prepare_data_for_images_example()