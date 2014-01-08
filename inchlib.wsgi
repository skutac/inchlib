#sys.stdout = sys.stderr
import os
import sys
import site

##Remember original sys.path.
#prev_sys_path = list(sys.path)
#
## Add each new site-packages directory.
#for directory in ALLDIRS:
#        site.addsitedir(directory)
#
## Reorder sys.path so new directories at the front.
#new_sys_path = []
#for item in list(sys.path):
#        if item not in prev_sys_path:
#                new_sys_path.append(item)
#                sys.path.remove(item)
#        print item
#sys.path[:0] = new_sys_path

os.environ['PYTHON_EGG_CACHE'] = '/tmp/python-eggs'

#If your project is not on your PYTHONPATH by default you can add the following
sys.path.append('/var/www/projects/inchlib/')
os.environ['DJANGO_SETTINGS_MODULE'] = 'inchlib.settings'

import django.core.handlers.wsgi
application = django.core.handlers.wsgi.WSGIHandler()
