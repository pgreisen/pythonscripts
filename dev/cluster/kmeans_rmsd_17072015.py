from pylab import plot,show
from numpy import vstack,array,median
from numpy.random import rand
from scipy.cluster.vq import kmeans,vq,kmeans2
from collections import defaultdict
import sys

def substract_mean(lst,key,centroids):
     min_value = 10000
     list_of_values = []
     for value in lst:
          #print type(value), type(centroids[key][0])
          print value - centroids[key][0]
          new_value = abs(value - centroids[key][0])
          if( new_value < min_value):
               min_value = value
     print "KEY is ",key
     print "CCCC",min_value
     print "DDDD",centroids[key][0]
     print "The closest point in cluster :",key
     print "RMSD = ", min_value
     if(min_value in lst):
          print min_value
          print "Found"

##########################
sort_cluster_points = defaultdict(list)
##########################

size_of_clusters = 8
tmpfile = sys.argv[1]
clusters = {}
data = ([])
with open(tmpfile,'r') as f:
     for line in f:
          tmpline = line.split()
          if(line[0] == "@" or line[0] == "&"):
               continue
          elif( len( tmpline ) == 2 ):
               clusters[ tmpline[1] ] = tmpline[0]
               tmp_rmsd = float( tmpline[1] )
               tmp_pair = array( [ tmp_rmsd, tmp_rmsd  ] )
               data.append( tmp_pair )
          else:
               continue

###############################################################
###############################################################
data = vstack(data)
centroids,labels = kmeans2( data, size_of_clusters )
#print labels
print centroids
#print "###################################"
#print len(data), len(labels), len(centroids)




for i in range( len(labels ) ):
     key = labels[i]
     sort_cluster_points[key].append(float(data[i][0]))
     

for key in sort_cluster_points:
     substract_mean(sort_cluster_points[key],key,centroids)

assert 1 ==0 

# assign each sample to a cluster
idx,_ = vq(data,centroids)
###help(kmeans)
###help(vq)
# some plotting using numpy's logical indexing
plot(data[idx==0,0],data[idx==0,1],'ob', data[idx==1,0],data[idx==1,1],'or')
plot(centroids[:,0],centroids[:,1],'sg',markersize=8)

show()

center_of_cluster = centroids[:,0]
for tm in center_of_cluster:
     print "###: ",tm
     if( tm in clusters.keys() ):
          print cluster[keys]
