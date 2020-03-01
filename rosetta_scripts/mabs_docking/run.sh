for i in $(seq 1 8);
	 do
./docking_protocol.static.linuxgccrelease -database database \
					  @flags -out:prefix $i -s $1 -native $1 > /dev/null & echo "done";
done
