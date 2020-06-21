hhblits -o /dev/null -i fasta/$dir/$id -oa3m hhblits/$dir/$id -d $db -n $iter -e $evalue -mact 0.35 -maxfilt 100000000 -neffmax 20 -cpu 16 -nodiff -realign_max 10000000 -maxmem 64
