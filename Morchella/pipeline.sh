#!/bin/bash
## genome assembly
# canu
## finisherSC
perl -pe 's/>[^\$]*$/">Seg" . ++$n ."\n"/ge' m54216_180521_075208.subreads.fasta > res/raw_reads.fasta
perl -pe 's/>[^\$]*$/">Seg" . ++$n ."\n"/ge' ydj.fasta > res/contigs.fasta
python /stor9000/apps/users/NWSUAF/2018055070/sf/finishingTool-2.1/finisherSC.py -par 20 res /stor9000/apps/users/NWSUAF/2018055070/sf/mummer/bin
## arrow
pbmm2 align -j 20 finisherSC/res/improved3.fasta bam.fofn | samtools sort > ydj_sorted.bam
pbindex ydj_sorted.bam
samtools faidx finisherSC/res/improved3.fasta
arrow ydj_sorted.bam -r finisherSC/res/improved3.fasta -o ydj_polished.fasta -j 20

## completeness assessment
/users/wanghm/miniconda2/bin/meryl threads=16 k=18 count m54216_180521_075208.subreads.fasta.gz output scls.meryl
/users/wanghm/sf/merqury/merqury.sh scls.meryl ../ydj.fasta scls_merqury
# coverage
minimap2 -ax map-pb -t 30 ydj.fasta m54216_180521_075208.subreads.fasta | samtools view -@ 10 -bS - | samtools sort -@ 10 -o ydj_pbSorted.bam -
samtools depth ydj_pbSorted_mapped.bam > ydj_pb.depth
python3 minBinDepth.py ydj_pb.depth ydj_pb.cov

## repeat annotation
BuildDatabase -name nztd_db -engine ncbi Morchella_NZTD.fasta
RepeatModeler -engine ncbi -pa 6 -database nztd_db
~/software/RepeatMasker/util/queryRepeatDatabase.pl -species fungi > fungi_repeats.lib
cat fungi_repeats.lib RM*/consensi.fa.classified > myGenome.custom.repeat.lib
RepeatMasker -engine ncbi -xsmall -nolow -lib myGenome.custom.repeat.lib -dir ./NZTD_ncbiSoft -gff -html -pa 4 -no_is Morchella_NZTD.fasta

## identity full-length LTR
GENOME=ydj.fasta
#LTRharvest
~/genometools-1.5.9/bin/gt suffixerator \
	  -db ${GENOME} \
	    -indexname ydj \
		  -tis -suf -lcp -des -ssp -sds -dna
~/genometools-1.5.9/bin/gt ltrharvest \
	  -index ydj \
	    -similar 85 -vic 10 -seed 20 -seqids yes \
		  -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 > ydj.harvest.scn
# LTR harvest Non-canonical LTR-RT candidates
~/genometools-1.5.9/bin/gt ltrharvest \
	 -index ydj  \
	    -similar 85 -vic 10 -seed 20 -seqids yes \
		 -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6  -motif TGCA -motifmis 1 > ydj.harvest.TGCA.scn

# LTR_FINDER
~/LTR_Finder/source/ltr_finder -D 15000 -d 1000 -L 7000 -l 100 -p 20 -C -M 0.85 ${GENOME} > ydj.finder.scn
# LTR_retriever integrate results of LTR_harvest and LTR_Finder
source activate LTR_retriever
~/LTR_retriever-master/LTR_retriever -genome ${GENOME} -inharvest ydj.harvest.TGCA.scn -infinder ydj.finder.scn -nonTGCA ydj.harvest.scn -threads 2

## gene model prediction
ref=ydj.fasta
MASKED_GENOME=ydj.fasta.masked
ortho_protein=/disk/alpha/mango/odb9_fungi/odb9_fungi.fasta
# 1. Augustus Abinitio prediction using braker2 integrating evidence from protein,RNA-seq
braker.pl --fungus --cores 12 --etpmode --softmasking --gff3 --genome=${MASKED_GENOME} --prot_seq=${ortho_protein} --bam=/disk/alpha/Morchella/gene_predict/ydj1.sorted.bam,/disk/alpha/Morchella/gene_predict/ydj2.sorted.bam
# 2. GeneMark_ES Abinitio prediction
perl ~/gm_et_linux_64/gmes_petap.pl --sequence ${ref} --ES --fungus --cores 14
# 3. GeneMark_ET integrate transcriptome
fq1=YDJ-1_HL7JHCCXY_L3_1.clean.fq.gz
fq2=YDJ-1_HL7JHCCXY_L3_2.clean.fq.gz
#mkdir -p star_index
#STAR \
#	    --runThreadN 14 \
#		    --runMode genomeGenerate \
#			    --genomeDir star_index \
#				    --genomeFastaFiles ${ref}
STAR \
    --runThreadN 6 \
		    --runMode alignReads \
			    --genomeDir star_index \
				     --readFilesIn ${fq1},${fq2} \
						--alignIntronMin 20 --alignIntronMax 4000 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outWigType wiggle read2
~/gm_et_linux_64/star_to_gff.pl --star  SJ.out.tab --gff STAR.SJ2.gff --label STAR

RNASEQ_hints=/disk/alpha/Morchella/gene_predict_STAR_hints/STAR.SJ2.gff
/home/wanghm/gm_et_linux_64/gmes_petap.pl --sequence ${ref} --ET ${RNASEQ_hints} --fungus --cores 8

# 4. extract candicate coding regions from denovo transcripts assembly
stringtie ${bam} --rf -p 14 -o ydj1.stingtie.gtf
trans_asm=ydj1.stringtie.gtf
~/EvidenceModeler/TransDecoder-TransDecoder-v5.5.0/util/gtf_genome_to_cdna_fasta.pl ${trans_asm} ${ref} > transcripts.fasta
/home/wanghm/EvidenceModeler/TransDecoder-TransDecoder-v5.5.0/util/gtf_to_alignment_gff3.pl ${trans_asm} > transcripts.gff3
/home/wanghm/EvidenceModeler/TransDecoder-TransDecoder-v5.5.0/TransDecoder.LongOrfs -t transcripts.fasta
/home/wanghm/EvidenceModeler/TransDecoder-TransDecoder-v5.5.0/TransDecoder.Predict -t transcripts.fasta
/home/wanghm/EvidenceModeler/TransDecoder-TransDecoder-v5.5.0/util/cdna_alignment_orf_to_genome_orf.pl transcripts.fasta.transdecoder.gff3 transcripts.gff3 transcripts.fasta > transcripts.fasta.transdecoder.genome.gff3

########### run Evmodeler using transDecoder output; braker2 output;geneMark_ES output as gene_predictions.gff3
##########                using geneMark _ET output as transcript_alignments.gff3

cat transcripts.fasta.transdecoder.genome.gff3 braker2.gff3 gene_predict_GeneMark_ES/genemark_ES.gff3 > gene_predictions.gff3
ln gene_predict_STAR_hints/genemark_ET.gff3 transcript_alignments.gff3
~/EvidenceModeler/EVidenceModeler-1.1.1/EvmUtils/partition_EVM_inputs.pl --genome ${ref} --gene_predictions gene_predictions.gff3 \
	--transcript_alignments transcript_alignments.gff3 --segmentSize 100000 --overlapSize 10000 --partition_listing partitions_list.out
~/EvidenceModeler/EVidenceModeler-1.1.1/EvmUtils/write_EVM_commands.pl --genome ${ref} --weights /disk/alpha/Morchella/evm_predict/weights.txt --gene_predictions gene_predictions.gff3 --transcript_alignments transcript_alignments.gff3 \
				      --output_file_name evm.out  --partitions partitions_list.out >  commands.list
~/EvidenceModeler/EVidenceModeler-1.1.1/EvmUtils/execute_EVM_commands.pl commands.list | tee run.log
~/EvidenceModeler/EVidenceModeler-1.1.1/EvmUtils/recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out
~/EvidenceModeler/EVidenceModeler-1.1.1/EvmUtils/convert_EVM_outputs_to_GFF3.pl  --partitions partitions_list.out --output evm.out  --genome ${ref}
find . -regex ".*evm.out.gff3" -exec cat {} \; > EVM.all.gff3
# reformat
gffread -T EVM.all.gff3 -o EVM.all.gtf
python3 ~/bin/gtfUtils.py -i EVM.all.gtf -reformat -o EVM.all.reformat.gtf
python3 ~/bin/gtfUtils.py -i EVM.all.reformat.gtf -r SCLS_ -o SCLS.gtf
gffread -y SCLS.protein.fasta -T SCLS.gtf -g OYR.fasta
awk '/>/{split($2,a,"=");print ">"a[2];next}{print}' SCLS.protein.fasta > SCLS_protein.fasta
rm SCLS.protein.fasta

## GO annotation
sh /stor9000/apps/users/NWSUAF/2018055070/database/interproscan-5.44-79.0/interproscan.sh -appl CDD,COILS,MobiDBLite,Gene3D,HAMAP,PANTHER,Pfam,PIRSF,PRINTS,ProDom,PROSITEPATTERNS,PROSITEPROFILES,SUPERFAMILY,TIGRFAM -i /stor9000/apps/users/NWSUAF/2018055070/data/Morchella/SCLS_protein.fasta -b /stor9000/apps/users/NWSUAF/2018055070/data/Morchella/scls_interpro -goterms -iprlookup -pa -f TSV -cpu 30 -dp
## SM clusters
python3 ~/PycharmProjects/test1/mango/antismash_gbk2fasta.py *.region*.gbk > scls_antiSMASHCoreGenes.fasta
awk '/>/{split($1,a,"[>;]");print a[2]}' scls_antiSMASHcoreGenes.fasta
for j in knownclusterblast/*.txt;do
	awk 'NR==FNR{if($3~/transcript/){FS="\"";id[$4]=$2}}NR>FNR{if($0 ~ /^evm/){FS=OFS="\t";print id[$1],$2,$3,$4}else{print}}' scls.gtf ${j} > ${j%.txt}_cluster.geneID
done
## CAZymes
awk '$6>=2{gsub("[0-9].*","",$2);gsub("\\.t.*","",$1);print $1,$2}' ydj_cazy.tbl | sort | uniq > ydj_cazy.geneID

# identify syntenic blocks
nucmer=/Users/alexwang/software/mummer4/bin/nucmer
${nucmer} --maxmatch -c 500 -b 500 -l 100 -t 8 --prefix=ydj nztd.fa ../ydj.fasta
delta-filter -m -i 95 -l 100 ydj.delta > ydj.delta.filter
show-coords -THrd ydj.delta.filter > ydj.coords
# LS regions
awk 'BEGIN{FS=OFS="\t"}{if($1<$2){print $10,$1,$2,"-",FILENAME}else{print $10,$2,$1,"-",FILENAME}}' ydj.coords | bedtools sort -i - | bedtools genomecov -i - -g nztd.genome.txt -bga > nztd.genomecov.txt
awk 'BEGIN{FS=OFS="\t"}{if($3<$4){print $11,$3,$4,"-",FILENAME}else{print $11,$4,$3,"-",FILENAME}}' ydj.coords | bedtools sort -i - | bedtools genomecov -i - -g ydj.genome.txt -bga > ydj.genomecov.txt
# snps
show-snps -ClrT ydj.delta.filter > nztd_ydj.snps

### find syntenic genes
awk 'BEGIN{FS=OFS="\t"}{$9="id="$9;print$0}' ../nztd.gtf > nztd_c.gtf
awk 'BEGIN{FS=OFS="\t"}{$9="id="$9;print$0}' ../SCLS.gtf > SCLS_c.gtf
python -m jcvi.formats.gff bed --type=gene --key=gene_id ../EVM.all.reformat.rename.gtf -o SCLS.bed
python -m jcvi.formats.gff bed --type=gene --key=id nztd_c.gtf -o nztd.bed
# this step must have cds fasta file in current directory, cds fasta file must have the same headers as bed file above
# find ortholog
python -m jcvi.compara.catalog ortholog scls nztd --cscore=.99
#build .simple
python -m jcvi.compara.synteny screen --minspan=10 --simple scls.nztd.anchors scls.nztd.anchors.new
awk -F ">" 'NR==FNR&&/>/{{printf$2","}}NR>FNR&&/>/{{printf$2","}}' ../ydj.fasta ../colin/nztd.fa > SCLS.nztd.seqid
## before plot, you should manual setting .seqid file to 2 lines
python -m jcvi.graphics.karyotype scls.nztd.seqid scls_nztd.layout telo.txt --figsize=5x5

## find orthogroups
python orthofinder.py -f OrthoFinder2 -M msa -T fasttree -t 14 -a 12
## select specie-specific gene families
awk '$2<$3 && $2!=0{print}' Orthogroups.GeneCount.tsv | cut -f 1 > scls.expand.orthogroups
awk '$2<$3 && $2==0{print}' Orthogroups.GeneCount.tsv | cut -f 1 > scls.sp.orthogroups
grep -f scls.expand.orthogroups Orthogroups.tsv | grep -o "SCLS....." scls.expand.geneID
grep -f scls.sp.orthogroups Orthogroups.tsv | grep -o "SCLS....." > scls.sp.geneID

## identify 6mA modification of fungi
pbmm2 align -j 30 ydj.fasta m54216_180521_075208.subreads.bam | samtools sort -o pbmm2_sorted.bam -
pbindex pbmm2_sorted.bam
/stor9000/apps/users/NWSUAF/2018055070/smrtlink/smrtcmds/bin/ipdSummary --reference ydj.fasta --gff ydj_mod_pbmm2.gff -j 16 --identify m6A,m4C --methylFraction --minCoverage 20 pbmm2_sorted.bam
ipdSummary --reference ydj.fasta --gff ydj_mod_pbmm2.gff -j 16 --identify m6A,m4C --methylFraction --minCoverage 20 pbmm2_sorted.bam
