import json, re

from pygments import highlight
from pygments.lexers import PythonLexer, JavascriptLexer, BashLexer
from pygments.formatters import HtmlFormatter

import config

from django.shortcuts import render_to_response, redirect
from django.utils.safestring import mark_safe

from examples.models import Examples, SettingsAttributes

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
    example.settings = highlight(settings, JavascriptLexer(), HtmlFormatter())
    return render_to_response("inchlib_examples.html", {"examples":examples, "example": example, "settings": settings, "next": next, "previous": previous})

def use_cases(req, exampleid):
    examples = [e for e in Examples.objects.filter(exampletype=3)]
    examples.sort(key=lambda e: e.order)
    example = Examples.objects.get(exampleid=exampleid)

    settings = example.examplesettings_set.all()
    settings = mark_safe(json.dumps(parse_settings({e.settingsattribute.name: e.value for e in settings})))

    example.description = re.sub('href="', '"'.join(['href=', config.BASE_URL]), example.description)
    return render_to_response("inchlib_use_cases.html", {"examples":examples, "example": example, "settings": settings})

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

def inchlib_clust(req):
    code = """
import inchlib_clust

c = inchlib_clust.Cluster()
c.read_csv(filename="filename", delimiter=",", header=bool)
# c.read_data(data, header=bool) use read_data() for list of lists instead of a data file
c.cluster_data(data_type="numeric", distance_measure="euclidean", linkage="ward", axis="both")

d = inchlib_clust.Dendrogram(c)
d.create_dendrogram(contract_clusters=bool, cluster_count=1000, write_data=bool)
d.add_metadata_from_file(metadata_file="filename", delimiter=",", header=bool)
d.export_dendrogram_as_json("filename")"""

    bash = "python inchlib_clust.py input_file.csv -m metadata.csv -dh -mh -d euclidean -l ward -a both -dd , -md ,"

    code = highlight(code, PythonLexer(), HtmlFormatter())
    bash = highlight(bash, BashLexer(), HtmlFormatter())

    return render_to_response("inchlib_clust.html", {"code": code, "bash": bash})
