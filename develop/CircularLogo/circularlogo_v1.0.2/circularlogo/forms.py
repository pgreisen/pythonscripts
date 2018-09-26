from django import forms
import json

ALPHABET_CHOICES = (
    ('DNA', ("DNA base")),
    ('RNA', ("RNA base"))
)

class MotifSiteForm(forms.Form):
    title = forms.CharField(required=False)
    inputText = forms.CharField(required=False, widget=forms.Textarea(attrs={'rows': 18, 'cols': 100}))
    inputFile = forms.FileField(required=False)    
    #alphabet = forms.ChoiceField(choices=ALPHABET_CHOICES,widget=forms.Select()) 
    height = forms.IntegerField(initial=500)
    width = forms.IntegerField(initial=500)
    pseudoA = forms.FloatField(required=False, initial=0.25, widget=forms.TextInput(attrs={'size': 8})) 
    pseudoT = forms.FloatField(required=False, initial=0.25, widget=forms.TextInput(attrs={'size': 8}))
    pseudoC = forms.FloatField(required=False, initial=0.25, widget=forms.TextInput(attrs={'size': 8}))
    pseudoG = forms.FloatField(required=False, initial=0.25, widget=forms.TextInput(attrs={'size': 8}))
    mCHOICES = (('mi', 'mutual information'), ('chi', 'Chi-square'),)
    method = forms.ChoiceField(required=False, choices=mCHOICES) 
    gbkgA   = forms.FloatField(required=False, initial=0.25, widget=forms.TextInput(attrs={'size': 8}))
    gbkgT = forms.FloatField(required=False, initial=0.25, widget=forms.TextInput(attrs={'size': 8}))
    gbkgC = forms.FloatField(required=False, initial=0.25, widget=forms.TextInput(attrs={'size': 8}))
    gbkgG = forms.FloatField(required=False, initial=0.25, widget=forms.TextInput(attrs={'size': 8}))
    pvalue  = forms.FloatField(required=False, initial=1e-5,min_value=0.0, max_value=1.00, widget=forms.TextInput(attrs={'size': 8}))
    
    fCHOICES=[('fasta','FASTA Sequences Format'), ('json','Customized JSON Motif Format')]
    format = forms.ChoiceField(choices=fCHOICES, widget=forms.RadioSelect(attrs={"onChange":"formatSelect();"}) )
    hidden_config = forms.CharField(required=False, initial="", widget=forms.HiddenInput())
    
    def clean_bak(self):
        if not (self.cleaned_data['inputFile'] or self.cleaned_data['inputText']):
            raise forms.ValidationError('Please input your sequences or json-format motif in text box or upload an appropriate file.')
        
        if(self.cleaned_data['inputText']):
            if(self.cleaned_data['format'] == 'fasta'):
                inputText = self.cleaned_data['inputText']  
                lines = inputText.strip().split('\n')
                error, lengthes, charset = False, [], set()
                for i in range(0, len(lines)):
                    if(i % 2 == 0):
                        if(not lines[i].startswith('>')):
                            error = True
                    else:
                        lengthes.append(len(lines[i].strip()))
                        for base in lines[i].strip():
                            charset.add(base)
                charset = sorted(list(charset))
                charset = ''.join(charset)  
                print("charset:", charset)                
                if(charset == 'ACGT'):
                    self.cleaned_data['ALPHABET'] = 'DNA'
                elif(charset == 'ACGU'):
                    self.cleaned_data['ALPHABET'] = 'RNA'
                else:
                     raise forms.ValidationError("Please make sure your input are correct RNA (ACGU) or DNA (ACGT) sequences!")
                if(len(lengthes) < 25):
                    raise forms.ValidationError("Please input enough sequences, >= 25!")            
                for k in lengthes[1:]:
                    if(k != lengthes[0]):
                        error = True
                if error:
                    raise forms.ValidationError("Please check the input sequences for correct FASTA format with the same length!")    
            else:  #validing the correct json format
                jsontext = self.cleaned_data['inputText']
                error, message = False, ''
                try:
                    parsed_json = json.loads(jsontext)
                except:
                    message = "Please provide the correct JSON format of motif!"
                    raise forms.ValidationError(message)  
                if(not 'nodes' in parsed_json):
                    error = True
                    message = "Please check the correct JSON format: no nodes in motif file!"
                else:
                    nodes = parsed_json['nodes']
                    if(len(nodes) == 0):
                        error = True
                        message = "Please check the correct JSON format: zero node in motif file!"
                    else:
                        for node in nodes:
                            for key in ['index', 'bit', 'base', 'freq', 'label']:
                                if( not key in node):
                                    error = True
                                    message = "Please check the correct JSON format: missed keys ('index', 'bit', 'base', 'freq', 'label') for node !"
                                    break    
                            if(node['freq']):
                                freq = node['freq']
                                if(sum(freq) > 1.0):
                                    error = True  
                                    message = "Please check the correct JSON format: sum of freq shoule be less than 1.0!"
                            if(error):
                                break      
                if(error):
                    raise forms.ValidationError(message)     
        #checking the uploaded input file            
        #elif(self.cleaned_data['inputFile']):
        #    print(self.cleaned_data['inputFile'])
        #    content = self.cleaned_data['inputFile'].read()
        #    print(content)
        return self.cleaned_data
        
    def clean(self):
        if(self.cleaned_data['inputText']):
            inputText = self.cleaned_data['inputText']  
        #checking the uploaded input file            
        elif(self.cleaned_data['inputFile']):
            #print(self.cleaned_data['inputFile'])
            inputText = self.cleaned_data['inputFile'].read()
        else:
            raise forms.ValidationError('Please input your sequences or json-format motif in text box or upload an appropriate file.')        
        if(self.cleaned_data['format'] == 'fasta'):
            inputText = inputText.decode("utf-8")
            lines = inputText.strip().split('\n')
            error, lengthes, charset = False, [], set()
            for i in range(0, len(lines)):
                if(i % 2 == 0):
                    if(not lines[i].startswith('>')):
                        error = True
                else:
                    lengthes.append(len(lines[i].strip()))
                    for base in lines[i].strip():
                        charset.add(base)
            charset = sorted(list(charset))
            charset = ''.join(charset)
            #print("charset:", charset)            
            if(charset == 'ACGT'):
                self.cleaned_data['ALPHABET'] = 'DNA'
            elif(charset == 'ACGU'):
                self.cleaned_data['ALPHABET'] = 'RNA'
            else:
                 raise forms.ValidationError("Please make sure your input are correct RNA (ACGU) or DNA (ACGT) sequences!")
            if(len(lengthes) < 20):
                raise forms.ValidationError("Please input enough sequences, >= 20!")            
            for k in lengthes[1:]:
                if(k != lengthes[0]):
                    error = True
            if error:
                raise forms.ValidationError("Please check the input sequences for correct FASTA format with the same length!")
        else:  #validing the correct json format
            error, message = False, ''
            try:
                parsed_json = json.loads(inputText)
            except:
                message = "Please provide the correct JSON format of motif!"
                raise forms.ValidationError(message)  
            if(not 'nodes' in parsed_json):
                error = True
                message = "Please check the correct JSON format: no nodes in motif file!"
            else:
                nodes = parsed_json['nodes']
                if(len(nodes) == 0):
                    error = True
                    message = "Please check the correct JSON format: zero node in motif file!"
                else:
                    for node in nodes:
                        for key in ['index', 'bit', 'base', 'freq', 'label']:
                            if( not key in node):
                                error = True
                                message = "Please check the correct JSON format: missed keys ('index', 'bit', 'base', 'freq', 'label') for node !"
                                break    
                        if(node['freq']):
                            freq = node['freq']
                            if(sum(freq) > 1.003):  # due to round issue, we loose a little bit of here
                                error = True  
                                message = "Please check the correct JSON format: sum of freq shoule be less than 1.0!"
                        if(error):
                            break      
            if(error):
                raise forms.ValidationError(message)   
                
        self.cleaned_data['inputText'] = inputText          
        return self.cleaned_data
       