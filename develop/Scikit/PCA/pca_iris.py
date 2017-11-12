from sklearn import decomposition
from sklearn import cluster, datasets
import pylab as pl
iris = datasets.load_iris()

pca = decomposition.PCA()

pca.fit(iris.data)
print(pca.explained_variance_)  

# help(pca)
pca = decomposition.PCA(n_components=2)

pca.fit(iris.data)

X = pca.transform(iris.data)

import pylab as pl
pl.scatter(X[:, 0], X[:, 1], c=iris.target) 
pl.show()
