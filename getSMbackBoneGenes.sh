python3 antismash_gbk2fasta.py *.region*.gbk > antismash_core_protein.fasta
awk '/>/{split($1,a,"[>;]");print a[2]}' antismash_core_protein.fasta > antismash_gene.id
