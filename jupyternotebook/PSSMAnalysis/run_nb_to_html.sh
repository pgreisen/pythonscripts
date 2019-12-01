echo "Parameter to be set is the date as output"
papermill /Users/pgreisen/pythonscripts/jupyternotebook/PSSMAnalysis/20180313_PSSM_analysis.ipynb pssm_analysis.ipynb -p date $1 -p STD_DEV $2
jupyter nbconvert --to html pssm_analysis.ipynb;
cp /Users/pgreisen/pythonscripts/jupyternotebook/PSSMAnalysis/Visualize_PSSM_w_matrix.ipynb .
jupyter nbconvert --to html Visualize_PSSM_w_matrix.ipynb --execute 