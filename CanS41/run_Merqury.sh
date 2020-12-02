#!/bin/bash
## run merqury
genome_size=59000000
./best_k.sh ${genome_size}
k=18
for i in fa1 fa2; do
    # 1. Build meryl dbs
    meryl threads=12 k=$k count ../subreads/${i}.gz output ${i}.meryl
done
# 2. Merge
meryl union-sum output CanS41.meryl fa*.meryl
~/sf/merqury/merqury.sh CanS41.meryl ../subreads/OYR.fasta CanS41_merqury

## run busco
for i in OYR 67-1 YKD0085; do
    python ~/busco/scripts/run_BUSCO.py -i /disk/alpha/OYR/${i}.fasta -o ${i}_busco_res \
         -m genome -l /disk/alpha/busco_dataset/fungi_odb9 -c 12
done




