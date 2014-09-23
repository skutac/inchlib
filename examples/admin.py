from django.contrib import admin
from examples.models import Examples, ExampleTypes, SettingsAttributes, ExampleSettings


class ExamplesAdmin(admin.ModelAdmin):
    list_display = ('title', 'exampletype', 'order',)
    fieldsets = (
        (None, {
            'classes': ('wide',),
            'fields': ('title', 'description', 'data', 'exampletype', 'order'),
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

class SettingsAttributesAdmin(admin.ModelAdmin):
    list_display = ('name', 'description', 'settingsattributetype')
    fieldsets = (
        (None, {
            'classes': ('wide',),
            'fields': ('name', 'description', 'settingsattributetype'),
        }),
    )
    class Media:
        css = {
            "all": ("css/admin.css",)
        }

admin.site.register(SettingsAttributes, SettingsAttributesAdmin)


class ExampleSettingsAdmin(admin.ModelAdmin):
    list_display = ('example', 'settingsattribute', 'value')
    fieldsets = (
        (None, {
            'classes': ('wide',),
            'fields': ('example', 'settingsattribute', 'value'),
        }),
    )
    class Media:
        css = {
            "all": ("css/admin.css",)
        }

admin.site.register(ExampleSettings, ExampleSettingsAdmin)