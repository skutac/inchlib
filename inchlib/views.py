import json, re, urllib, csv, os, copy

from pygments import highlight
from pygments.lexers import PythonLexer, JavascriptLexer, BashLexer
from pygments.formatters import HtmlFormatter

from django.shortcuts import render_to_response, redirect
from django.utils.safestring import mark_safe
from django.http import HttpResponse
from django.conf import settings

from examples.models import Examples, SettingsAttributes

from inchlib_forms import InteractiveExampleForm

import inchlib_clust as clust

try:
    with open(os.path.join(settings.ROOT, "static/source_data/proteins_report.csv"), "r") as pdb_input:
        reader = csv.DictReader(pdb_input, delimiter=",")
        PDB2DATA = {p["PDB ID"]:p for p in reader}

    with open(os.path.join(settings.ROOT, "static/source_data/proteins_go_names.csv"), "r") as go_names_input:
        reader = csv.DictReader(go_names_input, delimiter=",")
        GO2NAME = {g["go"]:g["name"] for g  in reader}

    with open(os.path.join(settings.ROOT, "static/source_data/example_data.csv"), "r") as input_data:
        example_data = [r for r in csv.reader(input_data, delimiter=",")]

    with open(os.path.join(settings.ROOT, "static/source_data/example_metadata.csv"), "r") as input_data:
        example_metadata = [r for r in csv.reader(input_data, delimiter=",")]

except Exception, e:
    print str(e)
    PDB2DATA = {}
    GO2NAME = {}

try:
    with open(os.path.join(settings.ROOT, "static/source_data/compound2scaffold.csv"), "r") as scaffolds:
        reader = csv.DictReader(scaffolds, delimiter=",")
        COMPOUND2SCAFFOLD = {p["chembl_id"]:p["scaffold"] for p in reader}

except Exception, e:
    print str(e)
    COMPOUND2SCAFFOLD = {}


def index(req):
    return render_to_response("inchlib_index.html", {})

def index_construction(req):
    return render_to_response("inchlib_construction.html", {})

def release_notes(req):
    return render_to_response("inchlib_release_notes.html", {})

def release_notes_inchlib_clust(req):
    return render_to_response("inchlib_release_notes_inchlib_clust.html", {})

def examples(req, exampleid):
    examples = [e for e in Examples.objects.filter(exampletype=1)]
    examples.sort(key=lambda e: e.order)
    example = Examples.objects.get(exampleid=exampleid)
    next, previous = get_neighbours(int(exampleid), examples)

    example_settings = example.examplesettings_set.filter()
    example_settings = mark_safe(json.dumps(parse_settings({e.settingsattribute.name: e.value for e in example_settings})))

    example.description = re.sub('href="', '"'.join(['href=', settings.BASE_URL]), example.description)

    template = "inchlib_examples.html"
    if exampleid == "18":
        template = "inchlib_examples_summary.html"
    elif exampleid == "5":
        template = "inchlib_examples_row_compression.html"

    return render_to_response(template, {"examples":examples, "example": example, "settings": example_settings, "next": next, "previous": previous})

def interactive_example(req):
    interactive_form = InteractiveExampleForm()
    return render_to_response("inchlib_interactive_example.html", {"interactive_form": interactive_form})

def use_cases(req, exampleid):
    examples = [e for e in Examples.objects.filter(exampletype=3)]
    examples.sort(key=lambda e: e.order)
    example = Examples.objects.get(exampleid=exampleid)
    
    example_settings = example.examplesettings_set.all()
    example_settings = mark_safe(json.dumps(parse_settings({e.settingsattribute.name: e.value for e in example_settings})))

    example.description = re.sub('href="', '"'.join(['href=', settings.BASE_URL]), example.description)
    template = "inchlib_use_cases.html"
    
    if exampleid == "16":
        template = "inchlib_use_cases_proteins.html"
    elif exampleid == "17":
        template = "inchlib_use_cases_whiskey.html"
    elif exampleid == "13":
        template = "inchlib_use_cases_chemical_biology.html"
    elif exampleid == "12":
        template = "inchlib_use_cases_microarrays.html"
    
    return render_to_response(template, {"examples":examples, "example": example, "settings": example_settings})

def docs(req):
    attributes = list(SettingsAttributes.objects.filter(settingsattributetype = 1))
    attributes.sort(key = lambda x: x.name)
    events = list(SettingsAttributes.objects.filter(settingsattributetype = 2))
    events.sort(key = lambda x: x.name)
    return render_to_response("inchlib_docs.html", {"attributes": attributes, "events":events})

def input_format(req):
    example = Examples.objects.get(exampletype=2)
    example.data = mark_safe(example.data)
    return render_to_response("inchlib_input_format.html", {"example": example})

def get_neighbours(exampleid, examples):
    ids = [e.exampleid for e in examples]
    current = ids.index(exampleid)

    next = False
    if current != len(ids) - 1:
        next = ids[current+1]
    
    previous = False
    if current > 0:
        previous= ids[current-1]
        
    return next, previous

def parse_settings(example_settings):
    for k in example_settings.keys():
        try:
            val = float(example_settings[k])
            example_settings[k] = val
        except Exception, e:
            if example_settings[k] == "True":
                example_settings[k] = True

            elif example_settings[k] == "False":
                example_settings[k] = False

            elif example_settings[k].startswith("["):
                example_settings[k] = list(example_settings[k].strip("[]").split(","))

    return example_settings


def download(req):
    return render_to_response("inchlib_download.html", {})

def contact(req):
    return render_to_response("inchlib_contact.html", {})

def performance(req):
    return render_to_response("inchlib_performance.html", {})

def dev(req):
    return render_to_response("inchlib_dev.html", {})

def inchlib_clust_doc(req):
    return render_to_response("inchlib_clust_doc.html", {})

def inchlib_examples_summary_html(req):
    return render_to_response("inchlib_examples_summary_html.html", {})

def inchlib_clust(req):
    code = """
import inchlib_clust

#instantiate the Cluster object
c = inchlib_clust.Cluster()

# read csv data file with specified delimiter, also specify whether there is a header row, the type of the data (numeric/binary) and the string representation of missing/unknown values
c.read_csv(filename="/path/to/file.csv", delimiter=",", header=bool, missing_value=str/False, datatype="numeric/binary")
# c.read_data(data, header=bool, missing_value=str/False, datatype="numeric/binary") use read_data() for list of lists instead of a data file

# normalize data to (0,1) scale, but after clustering write the original data to the heatmap
c.normalize_data(feature_range=(0,1), write_original=bool)

# cluster data according to the parameters
c.cluster_data(row_distance="euclidean", row_linkage="single", axis="row", column_distance="euclidean", column_linkage="ward")

# instantiate the Dendrogram class with the Cluster instance as an input
d = inchlib_clust.Dendrogram(c)

# create the cluster heatmap representation and define whether you want to compress the data by defining the maximum number of heatmap rows, the resulted value of compressed (merged) rows and whether you want to write the features
d.create_cluster_heatmap(compress=int, compressed_value="median", write_data=bool)

# read metadata file with specified delimiter, also specify whether there is a header row
d.add_metadata_from_file(metadata_file="/path/to/file.csv", delimiter=",", header=bool, metadata_compressed_value="frequency")

# read column metadata file with specified delimiter, also specify whether there is a 'header' column
d.add_column_metadata_from_file(column_metadata_file="/path/to/file.csv", delimiter=",", header=bool)

# export the cluster heatmap on the standard output or to the file if filename specified
d.export_cluster_heatmap_as_json("/home/ctibor/Desktop/to_delete.json")
#d.export_cluster_heatmap_as_html("/path/to/directory") function exports simple HTML page with embedded cluster heatmap and dependencies to given directory 
"""

    bash = "python inchlib_clust.py input_file.csv -m metadata.csv -cm column_metadata.csv -dh -mh -cmh -d euclidean -l ward -a both -dd , -md , -cmd ,"

    code = highlight(code, PythonLexer(), HtmlFormatter())
    bash = highlight(bash, BashLexer(), HtmlFormatter())

    return render_to_response("inchlib_clust.html", {"code": code, "bash": bash})

def get_pdb_file(req):
    keys = ["PDB ID","Chain ID","Structure Title","Resolution","Classification","Source","Biological Process","Cellular Component","Molecular Function","PubMed ID","Mesh Terms","DOI","Sequence","Chain Length"]
    pdb_id = req.GET["pdb_id"]
    webgl = True if str(req.GET["webgl"]) == "true" else False
    multiple = ["Biological Process", "Cellular Component", "Molecular Function"]
    pdb_data = {}

    if pdb_id in PDB2DATA:
        pdb_data = copy.deepcopy(PDB2DATA[pdb_id])

        for key in multiple:
            if pdb_data[key]:
                values = pdb_data[key].split("*")
                go = values[1].strip(" \n")
                pdb_data[key] = [go, GO2NAME[go]] 
            else:
                pdb_data[key] = ["", ""]

    if webgl:
        pdb_data["pdb_file"] = fetch_pdb(pdb_id)

    return HttpResponse(json.dumps(pdb_data))

def fetch_pdb(id):
  url = 'http://www.rcsb.org/pdb/files/%s.pdb' % id
  return urllib.urlopen(url).read()

def get_scaffolds(req):
    compounds = req.GET.getlist("compounds[]")

    scaffold2compound = {}
    for chembl_id in compounds:
        scaffold = COMPOUND2SCAFFOLD[chembl_id]
        if scaffold in scaffold2compound:
            scaffold2compound[scaffold].append(chembl_id)
        else:
            scaffold2compound[scaffold] = [chembl_id]

    scaffolds = [(s, scaffold2compound[s]) for s in scaffold2compound]
    scaffolds.sort(key = lambda x: len(x[1]), reverse=True)

    return HttpResponse(json.dumps(scaffolds))

def get_compressed_rows_json_by_node(req):
    row_ids = req.GET.getlist("row_ids[]")
    data = [example_data[0]]
    data.extend([r for r in example_data if r[0] in row_ids])
    print data

    c = clust.Cluster()
    c.read_data(data, header=True)
    c.cluster_data(data_type="numeric", row_distance="euclidean", row_linkage="ward", axis="both")

    d = clust.Dendrogram(c)
    d.create_cluster_heatmap(compress=20, compressed_value="median")
    d.add_metadata(example_metadata, header=True, metadata_compressed_value="median")
    json = d.export_cluster_heatmap_as_json()
    return HttpResponse(json)
