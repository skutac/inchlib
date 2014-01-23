from django.conf.urls.defaults import *
# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()

urlpatterns = patterns('',
    # Examples:
    url(r'^$', 'inchlib.views.index', name='index'),
    url(r'^examples/(\d+)$', 'inchlib.views.examples', name='examples'),
    url(r'^docs$', 'inchlib.views.docs', name='docs'),
    url(r'^input_format$', 'inchlib.views.input_format', name='input_format'),
    url(r'^download$', 'inchlib.views.download', name='download'),
    url(r'^contact$', 'inchlib.views.contact', name='contact'),
    # url(r'^inch/', include('inchlib.foo.urls')),

    # Uncomment the admin/doc line below to enable admin documentation:
    # url(r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    url(r'^admin/', include(admin.site.urls)),
)
