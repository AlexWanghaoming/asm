# find SNP
nucmer --mum --prefix=nucmer ../ydj.fasta ../Morchella_NZTD.fasta 
# keep best hits in reference but not in query
delta-filter -r nucmer.delta > nucmer.filter.delta
show-coords -d -l -r -L 10000 -I 95 -c -T nucmer.filter.delta > nucmer.coords
show-snps -Clr -T -x 1 nucmer.filter.delta | awk '$2~/[ATCG]/ && $3~/[ATCG]/{print$0}'> ydj.snp
#python3 MUMmerSNPs2VCF.py ydj.snp ydj.snp.vcf
