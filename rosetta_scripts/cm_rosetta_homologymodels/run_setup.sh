dst=20181223_rosetta_cm;
cd 20181223_rosetta_cm;
for i in $(seq 1 100);
do
    qsub -v path=$dst,prefix=$i qsub.sh;
done