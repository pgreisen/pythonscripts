# mutfile = "muts.txt"
# clusters = "cluster_centers.fasta"
papermill /Users/pgreisen/pythonscripts/jupyternotebook/AnalyseMutationsClusters/AnalyzeMutationalClusters.ipynb cluster_analysis.ipynb  -p mutfile $1 -p clusters $2
jupyter nbconvert --to html cluster_analysis.ipynb