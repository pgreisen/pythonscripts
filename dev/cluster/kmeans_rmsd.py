from pylab import plot,show
from numpy import vstack,array,median
from numpy.random import rand
from scipy.cluster.vq import kmeans,vq
import sys


tmpfile = sys.argv[1]

data = ([])


# 0-1
a = []
# 1-2
b = []
# 2-3
c = []
# 3-4
d = []
# 4-
e = []

with open(tmpfile,'r') as f:
     for line in f:
          if(line[0] == "@" or line[0] == "&"):
               continue
          elif(len(line.split()) == 2 ):
               tmp_rmsd = float( line.split()[1])
               tmp_pair = array( [ tmp_rmsd, tmp_rmsd  ] )
               data.append( tmp_pair )


               if( tmp_rmsd < 1.0 ):
                    a.append( tmp_rmsd )
               elif( tmp_rmsd < 2.0 ):
                    b.append(tmp_rmsd)
               elif(tmp_rmsd < 3.0):
                    c.append(tmp_rmsd)
               elif(tmp_rmsd < 4.0 ):
                    d.append( tmp_rmsd)
               else:
                    e.append( tmp_rmsd )

          else:
               continue

#X = []

# a = array(data)

# data generation
##data = vstack((rand(150,2) + array([.5,.5]),rand(150,2)))

data = vstack(data)

print a[len(a)/2]
print b[len(b)/2]
print c[len(c)/2]
print d[len(d)/2]
print e[len(e)/2]


# computing K-Means with K = 2 (2 clusters)
##centroids,_ = kmeans(data,2)
centroids,_ = kmeans(data,8)
# assign each sample to a cluster
##idx,_ = vq(data,centroids)

idx,_ = vq(data,centroids)

# some plotting using numpy's logical indexing
#plot(data[idx==0,0],data[idx==0,1],'ob', data[idx==1,0],data[idx==1,1],'or')

plot(data[idx==0,0],data[idx==0,1],'ob', data[idx==1,0],data[idx==1,1],'or')
plot(centroids[:,0],centroids[:,1],'sg',markersize=8)
##print centroids[:,0],centroids[:,1]
show()
