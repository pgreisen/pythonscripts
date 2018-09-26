from django.http import HttpResponse

from django.shortcuts import render
from circularlogo import forms
from circularlogo import models

def my_homepage_view(request):
    form = forms.MotifSiteForm( initial={'format':'fasta'} ) 
    if request.method == 'POST':
        form = forms.MotifSiteForm(request.POST, request.FILES)
        if(form.is_valid()):
            cd = form.cleaned_data
            ##print("Alphabet:", cd['ALPHABET']) 
            if(cd['format'] == 'fasta'):
                motifGraph, mL, mId = models.wrapFastaMotifModel(cd)
            else:
                motifGraph, mL, mId = models.wrapJsonMotifModel(cd)     
            visParam = wrapWebVisParam(cd, mL, mId)
            if(cd['format'] == 'fasta'):
                visParam['pvalue'] = cd['pvalue']
            #print(visParam)   
            #print(motifGraph)              
            return render(request, 'display.html', {'visParam': visParam, 'motifGraph': motifGraph})
    return render(request, 'index.html', {'form': form})  

def wrapWebVisParam(form_cd, mL, mId):
    visParam = {}
    visParam['height'] = form_cd['height']
    visParam['width']  = form_cd['width']
    visParam['title'] = str(mId)
    visParam['motifL'] = mL
    visParam['looptime'] = list(range(1, mL+1))
    visParam['hconfig'] = str(form_cd['hidden_config'])
    visParam['method'] = str(form_cd['method'])
    return visParam
    
        
            
#def my_helppage_view(request):
    #return render(request, 'help.html')      
#    return render(request, 'http://circularlogo.sourceforge.net/')     


     
    
