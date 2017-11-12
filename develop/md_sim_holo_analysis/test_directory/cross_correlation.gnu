# The script is executed like this
#  gnuplot -geometry 800x200 -persist cross_correlation.gnu
set output "cross_correlation.png";
set title "Cross Correlation";
set xlabel "Residue Index";
set ylabel "Residue Index";
# set term png small;
set terminal png
set view 0,0,1,1;
set pm3d map;
set cbrange [-1:1];
set palette defined (0 0 0 0, 1 0 0 1, 3 0 1 0, 4 1 0 0, 6 1 1 1);
splot "cross_correlation.dat" matrix;
#cmd = 'gnuplot -geometry 800x200 -persist '+ fname + '.gnuplot'
#    failure,output = commands.getstatusoutput( cmd )
