#!/bin/bash

for i in `seq $1 $2`;
do
    echo "Job $i"
    ./bin/main63.exe > out/main63_$i.log
    ./bin/parsePylist.exe out/main63_$i.log /data/cmcginn/GeneratorsHEPMC/PYTHIA6/5p02TeV/pthat500/outFile_$i.root
    rm out/main63_$i.log
done

echo "Complete"