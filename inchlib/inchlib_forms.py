
from django import forms


class InteractiveExampleForm(forms.Form):
    heatmap_colors = forms.CharField(label="Data color scale", max_length=10)
    max_quantile = forms.IntegerField(label="Max quantile value")
    min_quantile = forms.IntegerField(label="Min quantile value")
    middle_quantile = forms.IntegerField(label="Middle quantile value")