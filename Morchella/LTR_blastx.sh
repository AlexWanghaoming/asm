makeblastdb -in RepeatPeps.lib -dbtype prot
blastx -query LTR_internal.fasta -db RepeatPeps.lib -num_threads 10 -outfmt 6 -max_target_seqs 1 -out blastx.txt
