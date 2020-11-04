for j in Chain{A..Z};
do
    echo $j;
    cd $j;
    for i in *pdb;
    do
	grep '^ATOM' $i >> tmp.pdb;
    	mv tmp.pdb $i;
    done
    cd ..;
done
echo "########## DONE ##############";