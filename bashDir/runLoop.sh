#!/bin/bash

for i in `seq 0 100`;
do
    echo "Job $i"
    ./bin/main63.exe > out/main63.log
    ./bin/parsePylist.exe out/main63.log /data/cmcginn/GeneratorsHEPMC/PYTHIA6/pthat50/outFile_$i.root
    rm out/main63.log
done