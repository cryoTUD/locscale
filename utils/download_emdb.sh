#!/bin/bash
for emd in "$@"
do
	EMDB_DIR="EMD-"$emd
	echo "Downloading directory: "$EMDB_DIR
	wget --quiet -r -nH -nc --cut-dirs=5  ftp://ftp.ebi.ac.uk//pub/databases/emdb/structures/"$EMDB_DIR"/ &
	PID=$!
	i=1
	sp="/|\\-/|\\-"
	echo -n ' '
	while [ -d /proc/$PID ]
	do
		  printf "\b${sp:i++%${#sp}:1}"
		  sleep 1
	done
done







