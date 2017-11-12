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
          new_value = abs(value - centroids[key][0])
          if( new_value < min_value):
               min_value = value
     print "The closest point in cluster :",key
     print "RMSD = ", min_value
     # todo
     # insert assert statement that value is present
     #if(min_value in lst):
     #     #print min_value
     #     #print "Found"

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
<<<<<<< HEAD
import pdb;pdb.set_trace()
centroids,labels = kmeans2( data, size_of_clusters )
<<<<<<< HEAD
#print labels
# print centroids
#print "###################################"
#print len(data), len(labels), len(centroids)

=======
>>>>>>> 1e4bf29a8bcc54ce5b6534d3b167be9b38596cd4
=======
centroids,labels = kmeans2( data, size_of_clusters )
>>>>>>> f9b28dd6c3e7b734c74be14d38ab2178e029fcd2

#print "Labels are: ", labels


for i in range( len(labels ) ):
     key = labels[i]
     sort_cluster_points[key].append(float(data[i][0]))
     
for key in sort_cluster_points:
     #print "Input key is: ", key
     substract_mean(sort_cluster_points[key],key,centroids)

# assign each sample to a cluster
idx,_ = vq(data,centroids)
# some plotting using numpy's logical indexing
plot(data[idx==0,0],data[idx==0,1],'ob', data[idx==1,0],data[idx==1,1],'or')
plot(centroids[:,0],centroids[:,1],'sg',markersize=8)

#show()
#
'''
center_of_cluster = centroids[:,0]
for tm in center_of_cluster:
     #print "###: ",tm
     if( tm in clusters.keys() ):
          print cluster[keys]
'''
