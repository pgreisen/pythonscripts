echo "Parameter to be set is the date as output"
pth=/Users/pgreisen/pythonscripts/jupyternotebook/AnalysePerAA;
pdbfile=$1;
pdbfilename=${1%.pdb};
echo "File is : ", $pdbfilename;
papermill $pth/AnalysisPerAA.ipynb analysis.ipynb -p pdbfile $pdbfile
jupyter nbconvert --to html analysis.ipynb;
