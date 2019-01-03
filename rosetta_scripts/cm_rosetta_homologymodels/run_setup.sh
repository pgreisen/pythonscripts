dst=./;
for i in $(seq 1 100);
do
    qsub -v path=$dst,prefix=$i qsub.sh;
done