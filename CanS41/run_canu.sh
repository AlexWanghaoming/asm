#!/bin/bash
# Canu version 2.1
canu -correct -p OYR -d OYR_canu_assembly genomeSize=55m useGrid=false maxThreads=24 corMhapSensitivity=normal minReadLength=2000 minOverlapLength=500 -nanopore /users/wanghm/0data/OYR/fa1.gz /users/wanghm/0data/OYR/fa2.gz

canu -trim \
        -p OYR -d OYR_canu_assembly maxThreads=24 genomeSize=55m useGrid=false minReadLength=2000 minOverlapLength=500\
        -corrected -nanopore /users/wanghm/0data/OYR/OYR_canu_assembly/OYR.correctedReads.fasta.gz

# if coverage > 90 set correctedErrorRate=0.013
canu -assemble \
        -p OYR -d OYR_canu_assembly maxThreads=24 genomeSize=55m useGrid=false correctedErrorRate=0.035\
        -corrected -nanopore /users/wanghm/0data/OYR/OYR_canu_assembly/OYR.trimmedReads.fasta.gz
 
## purge haplotigs: use purged.fa for downstream analysis
asm=/users/wanghm/0data/OYR/OYR_canu_assembly/OYR.contigs.fasta
reads=/users/wanghm/0data/OYR/OYR_subreads.fa.gz

~/miniconda3/bin/minimap2 -t 20 -x map-ont ${asm} ${reads} | gzip -c - > aln.paf.gz
## stats paf
/users/wanghm/sf/purge_dups/bin/pbcstat aln.paf.gz
/users/wanghm/sf/purge_dups/bin/calcuts PB.stat > cutoffs 2> calcults.log
## Split an assembly if NNN in the assembly
/users/wanghm/sf/purge_dups/bin/split_fa $asm > asm.split
## do a self-self alignment
/users/wanghm/miniconda3/bin/minimap2 -t 20 -x asm5 -DP asm.split asm.split | gzip -c > asm.split.self.paf.gz
# purge haplotigs and overlap
/users/wanghm/sf/purge_dups/bin/purge_dups -2 -T cutoffs -c PB.base.cov asm.split.self.paf.gz > dups.bed 2> purge_dups.log
# get the purged primary and haplotigs sequences from draft assembly
/users/wanghm/sf/purge_dups/bin/get_seqs dups.bed $asm

## Polish asm using ont raw subreads
# with Racon software
genome=/users/wanghm/0data/OYR/subreads/OYR.genome.fa
# Iteration 1
/users/wanghm/miniconda3/bin/minimap2 -t 20 /users/wanghm/0data/OYR/subreads/OYR.genome.fa /users/wanghm/0data/OYR/subreads/OYR_subreads.fa.gz > ONTmin_IT0.paf
time /users/wanghm/miniconda3/bin/racon -t 20 /users/wanghm/0data/OYR/subreads/OYR_subreads.fa.gz ONTmin_IT0.paf /users/wanghm/0data/OYR/subreads/OYR.genome.fa > ONTmin_IT1.fasta
# Iteration 2
/users/wanghm/miniconda3/bin/minimap2 -t 20 ONTmin_IT1.fasta /users/wanghm/0data/OYR/subreads/OYR_subreads.fa.gz > ONTmin_IT1.paf
time /users/wanghm/miniconda3/bin/racon -t 20 /users/wanghm/0data/OYR/subreads/OYR_subreads.fa.gz ONTmin_IT1.paf ONTmin_IT1.fasta > ONTmin_IT2.fasta
# Iteration 3
/users/wanghm/miniconda3/bin/minimap2 -t 20 ONTmin_IT2.fasta /users/wanghm/0data/OYR/subreads/OYR_subreads.fa.gz > ONTmin_IT2.paf
time /users/wanghm/miniconda3/bin/racon -t 20 /users/wanghm/0data/OYR/subreads/OYR_subreads.fa.gz ONTmin_IT2.paf ONTmin_IT2.fasta > ONTmin_IT3.fasta
