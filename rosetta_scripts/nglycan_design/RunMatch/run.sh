for i in *.pos;
do
    dr=${i%.pos};
    mkdir $dr;
    cd $dr;
    cp ../files/* .;
    mv ../$i .;
    ln -s $i pos.pos;
    ../match.static.linuxgccrelease -database ~/database @matchflags;
    cd ..;
done
