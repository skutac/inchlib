from django.conf.urls.defaults import *
from django.views.generic.base import RedirectView

# Uncomment the next two lines to enable the admin:
from django.contrib import admin
from django.conf import settings
admin.autodiscover()

urlpatterns = patterns('',
    # # Examples:
    url(r'^$', RedirectView.as_view(url="/".join([settings.BASE_URL, 'home/'])), name="redirect"),
    url(r'^home/$', 'inchlib.views.index', name='index'),
    url(r'^examples/(\d+)$', 'inchlib.views.examples', name='examples'),
    url(r'^use_cases/(\d+)$', 'inchlib.views.use_cases', name='use_cases'),
    url(r'^docs$', 'inchlib.views.docs', name='docs'),
    url(r'^input_format$', 'inchlib.views.input_format', name='input_format'),
    url(r'^inchlib_clust$', 'inchlib.views.inchlib_clust', name='inchlib_clust'),
    url(r'^download$', 'inchlib.views.download', name='download'),
    url(r'^contact$', 'inchlib.views.contact', name='contact'),
    url(r'^performance$', 'inchlib.views.performance', name='performance'),
    url(r'^inchlib_clust_doc$', 'inchlib.views.inchlib_clust_doc', name='inchlib_clust_doc'),
    url(r'^inchlib_examples_summary_html$', 'inchlib.views.inchlib_examples_summary_html', name='inchlib_examples_summary_html'),
    url(r'^dev$', 'inchlib.views.dev', name='dev'),
    url(r'^mirna$', 'inchlib.views.mirna', name='mirna'),
    url(r'^affps$', 'inchlib.views.affps', name='affps'),
    url(r'^interactive_example$', 'inchlib.views.interactive_example', name='interactive_example'),
    url(r'^release_notes$', 'inchlib.views.release_notes', name='release_notes'),
    url(r'^release_notes_inchlib_clust$', 'inchlib.views.release_notes_inchlib_clust', name='release_notes_inchlib_clust'),
    url(r'^get_pdb_file$', 'inchlib.views.get_pdb_file', name='get_pdb_file'),
    url(r'^get_scaffolds$', 'inchlib.views.get_scaffolds', name='get_scaffolds'),
    url(r'^get_compressed_rows$', 'inchlib.views.get_compressed_rows_json_by_node', name='get_compressed_rows'),
    # url(r'^inch/', include('inchlib.foo.urls')),

    # Uncomment the admin/doc line below to enable admin documentation:
    # url(r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    url(r'^admin/', include(admin.site.urls)),
)
