import csv, urllib, re

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

def calculate_descriptors(pdbs):
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
        DesObject=PyPro.GetProDes(sequence) ##construct a GetProDes object
        descs = DesObject.GetAAComp() ##calculate 20 amino acid composition descriptors
        pdb2desc[p["PDB ID"]] = translate_descs(descs)
    return pdb2desc

def export_pdb2desc(pdb2desc):
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


def read_pdb_report(filename):
    with open(filename, "r") as input_file:
        pdbs = [r for r in csv.DictReader(input_file, delimiter=",")]

    return pdbs

def translate_descs(descs):
    return {char2aa[d]:descs[d] for d in descs}

def get_go_name(go):
    url = 'http://amigo.geneontology.org/amigo/term/{}'.format(go)
    page = urllib.urlopen(url).read()
    name = re.search("<h1>(.*?)</h1>", page).group(1)
    print go, name
    return name

def get_go_names(pdbs):
    go_keys = ["Biological Process", "Cellular Component", "Molecular Function"]
    go2name = {}

    for p in pdbs:
        for go in go_keys:
            values = p[go].split("*")
            if len(values) > 1:
                go_key = values[1].strip(" \n")
                if not go_key in go2name:
                    go2name[go_key] = get_go_name(go_key)
    return go2name

def export_go_names(go2name):    
    with open("../source_data/proteins_go_names.csv", "w") as output:
        output.write("go,name\n")

        for go in go2name:
            output.write(",".join([go, go2name[go]]))
            output.write("\n")
    return

if __name__ == '__main__':
    pdbs = read_pdb_report("../source_data/proteins_report.csv")
    go2name = get_go_names(pdbs)
    export_go_names(go2name)
    # pdb2desc = calculate_descriptors(pdbs)
    # export_pdb2desc(pdb2desc)

