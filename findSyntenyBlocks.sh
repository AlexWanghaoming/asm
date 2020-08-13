#!/bin/bash

nucmer=/Users/alexwang/software/mummer4/bin/nucmer
minMatchLength=1000
        echo "step1: run nucmer"
        ${nucmer} --maxmatch --maxgap=500 --mincluster=100 -t 6 --prefix=ydj nztd.fa ../ydj.fasta
        delta-filter -r ydj.delta > ydj.delta.filter
        echo "show-coords -d -l -L ${minMatchLength} -r -c -T ydj.delta.filter > ydj.tab.txt"
        show-coords -d -l -L ${minMatchLength} -r -c -T ydj.delta.filter > ydj.tab.txt
        
        echo "step2: bundle synteny blocks and assign telomere. "
        echo "python3 mummer_bundle_telo.py ydj.tab.txt ../ydj.telo.txt"
        python3 mummer_bundle_telo.py ydj.tab.txt ../ydj.telo.txt
