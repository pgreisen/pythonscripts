for j in {F..Z};
do 
    for i in *pdb; 
    do
	sh ~/bin/clean_rosetta_protein.sh $i $j; 
    done
    
    mkdir "Chain$j";
    mv *pdb$j.pdb;
done

