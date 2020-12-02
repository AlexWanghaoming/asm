#!/bin/bash
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

