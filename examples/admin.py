from django.contrib import admin
from examples.models import Examples, ExampleTypes


class ExamplesAdmin(admin.ModelAdmin):
    # fields = ('title', 'description')
    list_display = ('title', 'exampletype', 'order',)
    fieldsets = (
        (None, {
        	'classes': ('wide',),
            'fields': ('title', 'description', 'data'),
        }),
    )
    class Media:
        css = {
            "all": ("css/admin.css",)
        }

admin.site.register(Examples, ExamplesAdmin)

class ExampleTypesAdmin(admin.ModelAdmin):
    # fields = ('title', 'description')
    list_display = ('exampltypeid', 'exampletypename')
    class Media:
        css = {
            "all": ("css/admin.css",)
        }

admin.site.register(ExampleTypes, ExampleTypesAdmin)

