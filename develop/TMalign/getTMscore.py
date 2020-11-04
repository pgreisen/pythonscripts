import sys, shutil, os, subprocess, operator
from collections import defaultdict


'''

for w in sorted(d, key=d.get, reverse=True):
  print w, d[w]


'''

path = './'
score = defaultdict(list)
files = os.listdir(path)
threshold = 0.50

for fl in files:
    if fl.endswith('logfile'):
        pdbid = fl#.split('_')[0]
        ids = pdbid #[0:4]
        # read file and get score
        with open(fl,'r') as f:
            for line in f:
                if( len( line) > 18 and line[0:8] == "TM-score"):
                    score[ids].append( float( line.split()[1]) )


'''

sort by chain 1 or 2

'''



news = sorted(score.items(), key=lambda (k, v): v[1])

for key, value in news:
    if( value[0] > threshold ):
        print key, value
                        
