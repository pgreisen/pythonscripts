from Bio.Seq import Seq
from Bio import SeqIO
from Bio import motifs
import json

import numpy as np
from scipy import stats
#from scipy.cluster import hierarchy
import random

def wrapJsonMotifModel(form_cd):
    if(form_cd['inputText']):
        motifGraph = form_cd['inputText']
    #else:
    #    motifGraph = json.load(form_cd['inputFile'])
    #    motifGraph = json.dumps(motifGraph)  #dump json into string
    #mL = motifGraph.count("index")
    jsonGraph = json.loads(motifGraph)
    mL = len(jsonGraph['nodes'])
    mId = jsonGraph['id'] if 'id' in jsonGraph else ''
    return motifGraph, mL, mId

def wrapFastaMotifModel(form_cd):
    RNA_flag = False
    if(form_cd['ALPHABET'] == 'RNA'):
        RNA_flag = True
    instances = []
    if(form_cd['inputText']):
        lines = form_cd['inputText'].strip().split('\n')
        sites = [ item.strip().upper() for item in lines if not item.startswith('>') and len(item.strip()) > 0 ]
        for item in sites:
            if(RNA_flag):
                item = item.replace('U', 'T')
            instances.append(Seq(item))
    #else:
    #    iter = SeqIO.parse(form_cd['inputFile'], 'fasta') 
    #    for item in iter:
    #        instances.append(item.seq.upper())             
    m = motifs.create(instances)   
    mId = str(form_cd['title'])
    motifGraph = '{ "id": "' + mId + '",'  
    pvalue = form_cd['pvalue']
    pseudoA = form_cd['pseudoA']
    pseudoT = form_cd['pseudoT']
    pseudoC = form_cd['pseudoC']
    pseudoG = form_cd['pseudoG']
    pseudocounts={'A':pseudoA, 'C': pseudoC, 'G': pseudoG, 'T': pseudoT}
    gbkgA = form_cd['gbkgA']
    gbkgT = form_cd['gbkgT']
    gbkgC = form_cd['gbkgC']
    gbkgG = form_cd['gbkgG']
    background = {'A':gbkgA, 'C': gbkgC, 'G': gbkgG, 'T': gbkgT}
    motifGraph += '"background":{"key":["A","T","C","G"], '
    motifGraph += '"val":[' + str(gbkgA) + ',' + str(gbkgT) + ',' + str(gbkgC) + ',' + str(gbkgG) + ']},'
    motifGraph += '"pseudocounts":{"key":["A","T","C","G"], '
    motifGraph += '"val":[' + str(pseudoA) + ',' + str(pseudoT) + ',' + str(pseudoC) + ',' + str(pseudoG) + ']}, '
    nodes = processNodes(m, pseudocounts, background, RNA_flag)
    #node section
    motifGraph += '''
                    "nodes": [
                 ''' 
    for node in nodes:
        motifGraph += node + ','
    motifGraph = motifGraph[0:-1] + ']'
    #edge section
    method = form_cd['method']
    #print("Method:", method)
    edges = processEdges(m, pseudocounts, pvalue, method)
    if(len(edges) > 0):
        motifGraph += ', "links": ['
        for edge in edges:
            motifGraph += edge + ','
        motifGraph = motifGraph[0:-1] + ']'
    motifGraph += '}'
    return motifGraph, len(m), mId 

    
def processNodes(m, pseudocounts, background, RNA_flag):
    m.pseudocounts = pseudocounts
    m.background = background
    pwm = m.counts.normalize(pseudocounts=m.pseudocounts)
    #https://en.wikipedia.org/wiki/Sequence_logo
    s, n = len(m.alphabet.letters), len(m.instances)
    en = (1.0/(np.log(2))) * ((s - 1.0)/(2 * n))
    keys = list(pwm.keys())
    values = [pwm[k] for k in keys]
    pwm_np = np.asarray(values)
    pwm_t = pwm_np.transpose()
    H = np.nan_to_num(-pwm_t*np.log2(pwm_t))
    R = np.log2(s) - (H.sum(axis=1) + en)
    nodes = []
    for i in range(0, len(m)):
        indeces = np.argsort(pwm_t[i])
        prob = [ round(pwm_t[i][k], 3) for k in indeces ] 
        #base = [ keys[k] for k in indeces ]
        base = [ keys[k].replace('T', 'U') if RNA_flag else keys[k] for k in indeces ]
        node = '{"index": ' + str(i) + ', "label": "' + str(i+1) + '"'
        #node+= ', "bit": ' + str(R[i])
        node+= ', "bit": ' + str( max(R[i], 0.0) ) #to avoid negative values in "bit"
        node+= ', "base": ' + str(base) 
        node+= ', "freq": ' + str(prob) + ' }'
        nodes.append(node)
    return nodes


def processEdges(m, pseudocounts, pvalue, method):
    #simuScores = simulation(m, pseudocounts, W = 1000)
    simuMIDict = loadMISimuVals()
    #print("simuMIDict:", len(simuMIDict))
    rMI_u, rMI_s = simuMIDict[10000]
    N = len(m.instances)
    lookKey = int(10*divmod(N, 10)[0])
    if(lookKey in simuMIDict):
        rMI_u, rMI_s = simuMIDict[lookKey]
    k = len(m)
    interMat = np.zeros((k, k))
    if(method == 'chi'):
        calcDependence = calcDependence_chiSquare
    else:
        calcDependence = calcDependence_mutualInfo   
    edges = []
    for i in range(0, k):
        for j in range(i+1, k):
            #score, p = calcDependence_mutualInfo(m, i, j, pseudocounts)
            #score, p = calcDependence_chiSquare(m, i, j, pseudocounts) 
            score, p = calcDependence(m, i, j, pseudocounts) 
            interMat[i][j] = score
            if(method != 'chi'): ##updating p for MI method
                folder = abs(score - rMI_u)/rMI_s
                p = min(1.0, 1.0/folder/folder)
                #print(p)
            if(p <= pvalue ):
                edge = '{"source":' + str(i) + ', "target":' + str(j) + ', "value":' + str(score) + '}'
                edges.append(edge)     
    #print(interMat)
    #print(interMat.max(), interMat.min())
    return edges 
    
    
_BASE_indMap = {'A':0, 'C':1, 'G':2, 'T':3}    
def calcDependence_mutualInfo(m, i, j, pseudocounts):
    s = len(m.alphabet.letters)
    tbl = np.zeros((s, s)) 
    #count table of two cross-talk columns
    #      A     C     G     T
    #  A   11    12    13    14
    #  C   21    22    23    24
    #  G   31    32    33    34
    #  T   41    42    43    44
    for item in m.instances:
        iBase, jBase = item[i].upper(), item[j].upper()
        tbl[_BASE_indMap[iBase]][_BASE_indMap[jBase]] += 1 
    # to avoid zeros in calculation 
    pseudo = np.mat([pseudocounts[k] for k in ['A','T','C','G']])
    pseudoM = pseudo.transpose()*pseudo
    tbl +=  np.array(pseudoM)    
    tbl = tbl/tbl.sum() #to probability   
    PXs = tbl.sum(axis=1)
    PYs = tbl.sum(axis=0)
    Rij = 0
    for k in range(0, 4):
        for l in range(0, 4):
            pxy = tbl[k][l]
            px, py = PXs[k], PYs[l]
            if(pxy > 0):
                Rij += pxy*np.log2(pxy/px/py)  
    return Rij, 0
    
def calcDependence_chiSquare(m, i, j, pseudocounts):
    s = len(m.alphabet.letters)
    tbl = np.zeros((s, s)) 
    #count table of two cross-talk columns
    #      A     C     G     T
    #  A   11    12    13    14
    #  C   21    22    23    24
    #  G   31    32    33    34
    #  T   41    42    43    44
    for item in m.instances:
        iBase, jBase = item[i].upper(), item[j].upper()
        tbl[_BASE_indMap[iBase]][_BASE_indMap[jBase]] += 1 
    # to avoid zeros in calculation 
    pseudo = np.mat([pseudocounts[k] for k in ['A','C','G','T']])
    pseudoM = pseudo.transpose()*pseudo  
    tbl +=  np.array(pseudoM)    
    #tbl = tbl/tbl.sum() #to probability   
    values = [ v for v in tbl.flatten() if v > 0 ]
    chi, p = stats.chisquare(values)
    chi = 1.0*chi/len(values)
    return chi, p    
    
    
def simulation(m, pseudocounts, W = 5000):
    #pseudocounts={'A':0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}
    N = len(m.instances)
    L = len(m)
    scores, pvals = [], []
    for i in range(0, W):
        instances = []
        for j in range(0, N):
            chain = ''
            for k in range(0, L):
                chain += random.choice('ATCG')
            instances.append(Seq(chain))    
        rndM = motifs.create(instances)
        chi, p = calcDependence_chiSquare(rndM, 0, 1, pseudocounts)
        scores.append(chi)
        pvals.append(p)
        #Rij, p = calcDependence_mutualInfo(rndM, 0, 1, pseudocounts)
        #scores.append(Rij)
    return scores        
                    

def loadMISimuVals():
    handle = open('MI_simulation_output.txt')
    simuMIDict = {}
    for line in handle:
        if(line.startswith('#')): continue
        tokens = line.strip().split('\t')  
        simuMIDict[int(tokens[0])] = (float(tokens[1]), float(tokens[2])) 
    return simuMIDict
    
##50      0.147200970645  0.0684281837768                    
##200     0.0337882681921 0.0154617859455
