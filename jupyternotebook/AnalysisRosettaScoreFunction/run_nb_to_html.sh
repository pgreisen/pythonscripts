echo "Parameter to be set is the input score file from Rosetta"
pth=/Users/pgreisen/pythonscripts/jupyternotebook/AnalysisRosettaScoreFunction;
scorefile=$1;
scorefilename=${1%.sc};
echo "File is : ", $scorefilename;
papermill $pth/20190709_OA_biosensor_optimization.ipynb analysis.ipynb -p inputfile $scorefile
jupyter nbconvert --to html analysis.ipynb;
jupyter notebook analysis.ipynb;