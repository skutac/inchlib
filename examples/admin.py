from django.contrib import admin
from examples.models import Examples


class ExamplesAdmin(admin.ModelAdmin):
    # fields = ('title', 'description')
    list_display = ('title',)
    fieldsets = (
        (None, {
        	'classes': ('wide',),
            'fields': ('title', 'description'),
        }),
    )
    class Media:
        css = {
            "all": ("admin/css/admin.css",)
        }

admin.site.register(Examples, ExamplesAdmin)
