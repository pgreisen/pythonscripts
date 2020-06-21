# copy files to working directory
# tgz file from S3 name can be different 
for i in *tgz;
do
    tar zxf $i;
    name=${i%_pssm*};
done


for i in *pssm;
do
    if [ -d $i ]; then
	name=$i;
    fi
done

echo $name;
echo $STD_DEV;
ls -ltrls ;
cp $name/* .;
python PSSM_analysis.py $STD_DEV
Rscript Visualize_PSSM_w_matrix.r
zip $name\_viz\_results.zip *.png *csv *.txt
##aws s3 cp $name\_viz\_results.zip s3://tmpfilefasta/
