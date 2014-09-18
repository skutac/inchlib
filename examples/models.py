from django.db import models
from django.contrib import admin

# Create your models here.

class Examples(models.Model):
    exampleid = models.AutoField(primary_key=True)
    title = models.CharField(max_length=255)
    description = models.TextField()
    data = models.TextField()
    exampletype = models.ForeignKey("ExampleTypes")
    order = models.IntegerField()

    def __unicode__(self):
        return self.title

class SettingsAttributes(models.Model):
    settingsattributeid = models.AutoField(primary_key=True)
    name = models.CharField(max_length=255)
    description = models.TextField()
    settingsattributetype = models.IntegerField()

    def __unicode__(self):
        return self.name

class ExampleSettings(models.Model):
    examplesettingsid = models.AutoField(primary_key=True)
    example = models.ForeignKey("Examples")
    settingsattribute = models.ForeignKey("SettingsAttributes")
    value = models.CharField(max_length=255)

class ExampleTypes(models.Model):
    exampltypeid = models.AutoField(primary_key=True)
    exampletypename = models.CharField(max_length=255)

    def __unicode__(self):
        return self.exampletypename
