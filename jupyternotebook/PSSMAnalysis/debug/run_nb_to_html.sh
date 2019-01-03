echo "Parameter to be set is the date as output"
papermill /Users/pgreisen/pythonscripts/jupyternotebook/PSSMAnalysis/20180313_PSSM_analysis.ipynb pssm_analysis.ipynb -p date $1
#jupyter nbconvert --to html pssm_analysis.ipynb;
#jupyter nbconvert /Users/pgreisen/pythonscripts/jupyternotebook/PSSMAnalysis/Visualize_PSSM_w_matrix.ipynb