# Negative results - should not produce any output
python get_cststring_and_pdbfile.py -f 4981_299_1tdu_design_.pdb -a PHE,TRP
echo "No results"
python get_cststring_and_pdbfile.py -f 4981_299_1tdu_design_.pdb -a ASP,TRP
echo "Should produce output"
python get_cststring_and_pdbfile.py -f 4981_299_1tdu_design_.pdb -a GLU,TRP
echo "Should not produce output"
python get_cststring_and_pdbfile.py -f 4981_299_1tdu_design_.pdb -a ASP,TRP -c True
echo "Should copy file to directory - substringmatch/"
