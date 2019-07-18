echo "Parameter to be set is the input score file from Rosetta"
pth=/Users/pgreisen/pythonscripts/jupyternotebook/AnalyzeDockingResultsRosetta;
nameofplot=$1;
echo "File is : ", $nameofplot;
papermill $pth/analyzedockingrosetta.ipynb analysis.ipynb -p var1 $nameofplot
jupyter nbconvert --to html analysis.ipynb;
# jupyter notebook analysis.ipynb;