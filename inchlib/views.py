import json

from django.shortcuts import render_to_response
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
	#settings = mark_safe(json.dumps(parse_settings({e.settingsattribute.name: e.value for e in settings})))

	data = example.data
	example.data = mark_safe(example.data)

	return render_to_response("inchlib_examples.html", {"examples":examples, "example": example, "settings": settings, "data":data, "next": next, "previous": previous})

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

