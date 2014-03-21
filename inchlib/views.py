import json, re, urllib, csv

from pygments import highlight
from pygments.lexers import PythonLexer, JavascriptLexer, BashLexer
from pygments.formatters import HtmlFormatter

import config

from django.shortcuts import render_to_response, redirect
from django.utils.safestring import mark_safe
from django.http import HttpResponse

from examples.models import Examples, SettingsAttributes

try:
    with open("inchlib/static/source_data/proteins_report.csv", "r") as pdb_input:
        reader = csv.DictReader(pdb_input, delimiter=",")
        PDB2DATA = {p["PDB ID"]:p for p in reader}
except Exception, e:
    PDB2DATA = {}


def index(req):
    return render_to_response("inchlib_index.html", {})

def examples(req, exampleid):
    examples = [e for e in Examples.objects.filter(exampletype=1)]
    examples.sort(key=lambda e: e.order)
    example = Examples.objects.get(exampleid=exampleid)
    next, previous = get_neighbours(int(exampleid), examples)

    settings = example.examplesettings_set.all()
    settings = mark_safe(json.dumps(parse_settings({e.settingsattribute.name: e.value for e in settings})))

    example.description = re.sub('href="', '"'.join(['href=', config.BASE_URL]), example.description)
    return render_to_response("inchlib_examples.html", {"examples":examples, "example": example, "settings": settings, "next": next, "previous": previous})

def use_cases(req, exampleid):
    examples = [e for e in Examples.objects.filter(exampletype=3)]
    examples.sort(key=lambda e: e.order)
    example = Examples.objects.get(exampleid=exampleid)
    
    settings = example.examplesettings_set.all()
    settings = mark_safe(json.dumps(parse_settings({e.settingsattribute.name: e.value for e in settings})))

    example.description = re.sub('href="', '"'.join(['href=', config.BASE_URL]), example.description)
    template = "inchlib_use_cases.html"
    
    if exampleid == "16":
        template = "inchlib_use_cases_proteins.html"
    elif exampleid == "17":
        template = "inchlib_use_cases_whiskey.html"
    elif exampleid == "13":
        template = "inchlib_use_cases_chemical_biology.html"
    
    return render_to_response(template, {"examples":examples, "example": example, "settings": settings})

def docs(req):
    attributes = list(SettingsAttributes.objects.all())
    attributes.sort(key = lambda x: x.name)
    return render_to_response("inchlib_docs.html", {"attributes": attributes})

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

def parse_settings(settings):
    for k in settings.keys():
        try:
            val = float(settings[k])
            settings[k] = val
        except Exception, e:
            if settings[k] == "True":
                settings[k] = True

            elif settings[k] == "False":
                settings[k] = False

            elif settings[k].startswith("["):
                settings[k] = list(settings[k].strip("[]").split(","))

    return settings


def download(req):
    return render_to_response("inchlib_download.html", {})

def contact(req):
    return render_to_response("inchlib_contact.html", {})

def performance(req):
    return render_to_response("inchlib_performance.html", {})

def test(req):
    return render_to_response("inchlib_test.html", {})

def inchlib_clust_doc(req):
    return render_to_response("inchlib_clust_doc.html", {})

def inchlib_clust(req):
    code = """
import inchlib_clust

#instantiate the Cluster object
c = inchlib_clust.Cluster()

# read csv data file with specified delimiter, also specify whether there is a header row
c.read_csv(filename="filename", delimiter=",", header=bool)
# c.read_data(data, header=bool) use read_data() for list of lists instead of a data file

# normalize data to (0,1) scale, but after clustering write the original data to the heatmap
c.normalize_data(feature_range=(0,1), write_original=True)

# cluster data according to the parameters
c.cluster_data(data_type="numeric", distance_measure="euclidean", linkage="ward", axis="both")

# instantiate the Dendrogram class with the Cluster instance as an input
d = inchlib_clust.Dendrogram(c)

# create the cluster heatmap representation and define whether you want to contract the data, how much and if you want to write the features
d.create_dendrogram(contract_clusters=bool, cluster_count=1000, write_data=bool)

# read metadata file with specified delimiter, also specify whether there is a header row
d.add_metadata_from_file(metadata_file="filename", delimiter=",", header=bool)

# export the dendrogram on the standard output or to the file if filename specified
d.export_dendrogram_as_json("filename")"""

    bash = "python inchlib_clust.py input_file.csv -m metadata.csv -dh -mh -d euclidean -l ward -a both -dd , -md ,"

    code = highlight(code, PythonLexer(), HtmlFormatter())
    bash = highlight(bash, BashLexer(), HtmlFormatter())

    return render_to_response("inchlib_clust.html", {"code": code, "bash": bash})

def get_pdb_file(req):
    keys = ["PDB ID","Chain ID","Structure Title","Resolution","Classification","Source","Biological Process","Cellular Component","Molecular Function","PubMed ID","Mesh Terms","DOI","Sequence","Chain Length"]
    pdb_id = req.GET["pdb_id"]
    pdb_data = {}
    if pdb_id in PDB2DATA:
        pdb_data = PDB2DATA[pdb_id]
        
    pdb_data["pdb_file"] = fetch_pdb(pdb_id)
    return HttpResponse(json.dumps(pdb_data))

def fetch_pdb(id):
  url = 'http://www.rcsb.org/pdb/files/%s.pdb' % id
  return urllib.urlopen(url).read()

# def parse_pdb(pdb_file):
#     """
#     HEADER    TRANSFERASE/TRANSFERASE INHIBITOR       14-SEP-11   3TTJ              
#     TITLE     CRYSTAL STRUCTURE OF JNK3 COMPLEXED WITH CC-359"," A JNK INHIBITOR FOR  
#     TITLE    2 THE PREVENTION OF ISCHEMIA-REPERFUSION INJURY 
#     REMARK   2 RESOLUTION.    2.10 ANGSTROMS.
#     SOURCE   2 ORGANISM_SCIENTIFIC: HOMO SAPIENS;
#     JRNL        PMID   22226655 
#     JRNL        DOI    10.1016/J.BMCL.2011.12.028
#     """
#     data = {}
#     classification = re.search("HEADER\s+(.*?)\d+", pdb_file)
#     name = re.findall("TITLE\s+(\d+\s+)?(.*)", pdb_file)
#     name = " ".join([n[1].strip(" ") for n in name])
#     return data
