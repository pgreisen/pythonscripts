path=/z/insulin/users/pjug/pythonscripts/dev/TMalign;
for i in *logfile; 
do
	python $path/getTMscore.py $i >> TMscores.dat;
done
