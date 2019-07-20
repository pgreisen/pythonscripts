echo "Parameter to be set is the input score file from Rosetta, score term1, score term2, variants to analyse"
pth=/Users/pgreisen/pythonscripts/jupyternotebook/AnalysisRosettaScoreFunction;
scorefile=$1;
scorefilename=${1%.sc};
xvalue=$2;
yvalue=$3;
var1=$4;
echo "File is : ", $scorefile;
echo "Score terms set as : ", $x, $y
papermill $pth/AnalyzeRosettaScorefile.ipynb analysis.ipynb -p inputfile $scorefile -p xvalue $xvalue -p yvalue $yvalue -p var1 $var1
jupyter nbconvert --to html analysis.ipynb;
jupyter notebook analysis.ipynb;