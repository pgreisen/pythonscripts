for i in *.pos;
do
    dr=${i%.pos};
    mkdir $dr;
    cd $dr;
    cp ../files/* .;
    mv ../$i .;
    ln -s pos.pos $i;
    ../match.static.linuxgccrelease -database ~/database @matchflags;
    cd ..;
done
