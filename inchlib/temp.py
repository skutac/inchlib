# -*- coding: utf-8 -*-
import csv, shutil, MySQLdb, re, os

import rdkit
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit import RDConfig

class FragmentAnalyzer():

    def __init__(self, smiles, substructure_smiles):
        self.mol = AllChem.MolFromSmiles(smiles)
        self.substructure_smiles = substructure_smiles
        self.substructure = AllChem.MolFromSmiles(substructure_smiles)
        # self.mol_to_skeleton_map = self.get_mol_to_skeleton_mapping()

    # def get_mol_to_skeleton_mapping(self):
    #     atom_map = re.findall(":(\d+)\]", self.mapped_skeleton_smi)
    #     match = self.mapped_skeleton.GetSubstructMatch(self.mapped_skeleton)
    #     mol_to_skeleton_map = dict([(match[x], int(atom_map[x])) for x in range(len(atom_map))])
    #     return mol_to_skeleton_map

    def get_fragments(self):
        fragments = Chem.ReplaceCore(self.mol, self.substructure, labelByIndex=True)
        if fragments is None:
            return False

        # fragments = self.translate_fragment_positions(fragments)
        sidechains = Chem.GetMolFrags(fragments, asMols=True)
        sidechains = [Chem.MolToSmiles(s, True) for s in sidechains]
        return sidechains

    def get_features(self):
        fragments = Chem.ReplaceCore(self.mol, self.substructure, labelByIndex=True)
        fragments_smiles = Chem.MolToSmiles(fragments, True)
        features = self.get_position_features()
        return features

    def get_position_features(self):
        features = []
        fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
        factory = Chem.ChemicalFeatures.BuildFeatureFactory(fdefName)
        mol_features = factory.GetFeaturesForMol(self.mol)
        for a in mol_features:
            for pos in a.GetAtomIds():
                # if pos in self.mol_to_skeleton_map.keys():
                  feature = (pos, a.GetFamily(), a.GetType())
                  features.append(feature)
                # else:
                #     print pos, a.GetFamily(), a.GetType()

        return features

    # def translate_fragment_positions(self, fragments):
    #     sidechains = Chem.GetMolFrags(fragments, asMols=True)
    #     fragments = []
    #     for s in sidechains:
    #         fragment = Chem.MolToSmiles(s, True)
    #         if re.match("\[.*?\*\]", fragment):
    #             position = fragment[0:fragment.index("*")].strip("[")
    #             if position:
    #                 position = int(position)
    #             else:
    #                 self.mol_to_skeleton_map[position] = self.mol_to_skeleton_map[0]

    #             position = self.mol_to_skeleton_map[position]
    #             fragment = fragment[fragment.index("]")+1:]
            
    #             if re.search("\[\d+\*\]", fragment):
    #                 connection = re.search("\[(\d+)\*\]", fragment).group(1)
    #                 fragment = re.sub(connection, str(self.mol_to_skeleton_map[int(connection)]), fragment)

    #             fragments.append((position, fragment))
                    
    #     return fragments


def get_cursor(db):
    conn = MySQLdb.connect(host = "localhost", user = "root", passwd = "gugoun", db = db)
    cursor = conn.cursor(cursorclass=MySQLdb.cursors.DictCursor)
    return cursor

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
    with open("static/source_data/steroid_fragments.csv", "r") as input_file:
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

def create_era_files(ligands):
    with open("static/source_data/era_fragments.csv", "r") as input_file:
        rows = [r for r in csv.reader(input_file)]

    chemblid2fps = {}
    header = rows[0]
    chemblid2fragments = {r[0]:r[1:] for r in rows[1:]}

    for chemblid, fs in chemblid2fragments.items():
        chemblid2fps[chemblid] = []
        for f in fs:
            if f != "":
                f = re.sub("hash", "#", f)
                chemblid2fps[chemblid].append(get_morgan_fingerprint_for_smiles(f, radius=2, length=256))
            else:
                chemblid2fps[chemblid].append("".join(["0" for i in range(256)]))
    

    with open("static/source_data/era_fragment_fps.csv", "w") as w:
        writer = csv.writer(w)
        unique = set()

        for chemblid, fp in chemblid2fps.items():
            fp = "".join(fp)
            if not fp in unique:
                unique.add(fp)
                row = [chemblid]
                row.extend([i for i in fp])
                writer.writerow(row)

    with open("static/source_data/era_pic50.csv", "w") as w:
        unique = set()
        writer = csv.writer(w)
        writer.writerow(["id", "pIC50"])
        for l in ligands:
            if l["chemblid"] in chemblid2fps and not l["chemblid"] in unique:
                writer.writerow([l["chemblid"], round(float(l["value"]), 3)])
                unique.add(l["chemblid"])

def get_ligands_fragments(ligands):
    scaffold = "c1ccc2c(c1)CCC1C3CCCC3CCC21"
    ligand2fragments = {}
    positions = set()
    for l in ligands:
        try:
            fa = FragmentAnalyzer(l["smiles"], scaffold)
            fragments = fa.get_fragments()

            if fragments:
                ligand2fragments[l["chemblid"]] = {}

                for f in fragments:
                    if re.match("\[.*?\*\]", f):
                        original_position = f[0:f.index("*")].strip("[")
                        position = int(original_position) + 1 if original_position != "" else 1
                        positions.add(position)
                        f = re.sub("{}\*".format(original_position), "{}*".format(position), f)
                        generate_mol_svg(f)
                        f = re.sub("#", "hash", f)
                        ligand2fragments[l["chemblid"]][position] = f

                        generate_mol_svg(l["smiles"], name="{}_subs".format(l["chemblid"]), path="era_molecules", ext="png", size=(300,300), substructure=scaffold)

        except Exception as e:
            print(e)

    data = []
    positions = list(positions)
    positions.sort()

    for l, f in ligand2fragments.items():
        row = [l]
        for p in positions:
            row.append(f.get(p, None))
        data.append(row)

    with open("static/source_data/era_fragments.csv", "w") as w:
        writer = csv.writer(w)
        header = ["id"]
        header.extend(["Position {}".format(p) for p in positions])
        writer.writerow(header)
        writer.writerows(data)

    return data    

def generate_mol_svg(smiles, name=False, path="fragments", ext="svg", size=(100, 50), substructure=False):
    
    molpath = os.path.join("static", "img", path)

    if not os.path.exists(molpath):
        os.mkdir(molpath)
    
    smiles_path = re.sub("/", "", smiles)
    smiles_path = re.sub("#", "hash", smiles_path)
    filepath = os.path.join(molpath, "{}.{}".format(name if name else smiles_path, ext))
    try:
        mol = Chem.MolFromSmiles(smiles)
        
        if substructure:
            highlight_obj = Chem.MolFromSmiles(substructure)
            AllChem.Compute2DCoords(highlight_obj)
            AllChem.GenerateDepictionMatching2DStructure(mol, highlight_obj)
            matching = mol.GetSubstructMatch(highlight_obj)
            Draw.MolToFile(mol, filepath, size=size, highlightAtoms=matching)
        else:
            Draw.MolToFile(mol, filepath, size=size)
    except Exception as e:
        print(smiles, e)

def get_target_ligands(tid=19):
    """
    Get activities for given target
    """
    print("""Target ID:""", tid)
    cursor = get_cursor("chembl_19")
    
    cursor.execute("""SELECT activities.activity_id,
        activities.assay_id,
        activities.molregno,
        activities.standard_relation ,
        activities.standard_value,
        activities.standard_units,
        activities.standard_flag,
        activities.standard_type,
        activities.activity_comment,
        activities.pchembl_value,
        assays.confidence_score,
        assays.assay_type,
        target_dictionary.target_type as tgt_type,
        target_dictionary.pref_name as tgt_pref_name,
        target_dictionary.chembl_id as tgt_chembl_id,
        target_dictionary.organism,
        molecule_dictionary.pref_name as cmpd_pref_name,
        compound_structures.canonical_smiles,
        molecule_dictionary.chembl_id AS cmpd_chembl_id,
        assays.tid
        FROM activities,
        assays,
        target_dictionary,
        target_components,
        molecule_dictionary,
        compound_properties,
        compound_structures
        WHERE activities.potential_duplicate is null
        AND activities.data_validity_comment is null
        AND activities.assay_id = assays.assay_id
        AND assays.tid = target_dictionary.tid
        AND target_dictionary.tid = target_components.tid
        AND molecule_dictionary.molregno = activities.molregno
        AND molecule_dictionary.molregno = compound_properties.molregno
        AND molecule_dictionary.molregno = compound_structures.molregno
        AND activities.standard_relation = '='
        AND activities.standard_type IN ('IC50', 'pIC50', 'log IC50', 'Log IC50', 'logIC50')
        AND assays.tid = {}""".format(tid))

    compounds = cursor.fetchall()
    print("Raw data:", len(compounds), "bioactives")
    return compounds

if __name__ == '__main__':
    # prepare_data_for_images_example()
    # create_fragment_images()
    scaffold = "c1ccc2c(c1)CCC1C3CCCC3CCC21"
    generate_mol_svg(scaffold, "steran", ext="png", size=(500, 300), path="")
    # ligands = [{"chemblid": l["cmpd_chembl_id"], "value": l["pchembl_value"], "smiles": l["canonical_smiles"]} for l in get_target_ligands()]
    # data = get_ligands_fragments(ligands)
    # create_era_files(ligands)
    # create_fragment_fps()
    # print get_morgan_fingerprint_for_smiles("[12*]O", radius=2, length=256)