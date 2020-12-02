#!/bin/bash
ref=/disk/alpha/OYR/OYR.fasta
#mkdir RNAseq/star_index
#STAR --runThreadN 14 --runMode genomeGenerate --genomeDir RNAseq/star_index --genomeFastaFiles ${ref}
array=(67-1_24a 67-1_48a 67-1_8a C-24-1 C-8-1 S-24-1 S-8-1 MeOH_rep1 DON_rep2 PG_rep1)
length=${#array[@]}
#for ((i=0;i<${length};i++))
#do
#	hisat2 -x RNAseq/index/OYR -p 14 --min-intronlen 20 --max-intronlen 4000 --no-unal --no-temp-splicesite \
#		--novel-splicesite-outfile RNAseq/${array[$i]}.splicesite --novel-splicesite-infile RNAseq/${array[$i]}.splicesite \
#		-1 RNAseq/${array[$i]}_1.fastq -2 RNAseq/${array[$i]}_2.fastq | samtools view -bS - | samtools sort - RNAseq/${array[$i]}.sorted
#	STAR \
#		--runThreadN 8 \
#		    --runMode alignReads \
#			    --genomeDir RNAseq/star_index \
#				     --readFilesIn RNAseq/${array[$i]}_1.fastq,RNAseq/${array[$i]}_2.fastq \
#						--alignIntronMin 20 --alignIntronMax 4000 --outFilterIntronMotifs RemoveNoncanonical \
#						--outSAMtype BAM SortedByCoordinate --outFileNamePrefix RNAseq/${array[$i]}
## denovo assemble transcript
#	stringtie RNAseq/${array[$i]}.sorted.bam -p 14 -o ${array[$i]}.stingtie.assembly.gtf
#done

MASKED_GENOME=/disk/alpha/OYR/OYR.fasta.masked
ortho_protein=/disk/alpha/mango/odb9_fungi/odb9_fungi.fasta
# 1. Augustus Abinitio prediction using braker2 integrating evidence from protein,RNA-seq
braker.pl --fungus --cores 14 --etpmode --softmasking --gff3 --genome=${MASKED_GENOME} --prot_seq=${ortho_protein} --bam=RNAseq/67-1_24a.sorted.bam,RNAseq/67-1_48a.sorted.bam,RNAseq/67-1_8a.sorted.bam,RNAseq/C-24-1.sorted.bam,RNAseq/S-24-1.sorted.bam,RNAseq/C-8-1.sorted.bam,RNAseq/S-8-1.sorted.bam,RNAseq/DON_rep2.sorted.bam,RNAseq/MeOH_rep1.sorted.bam,RNAseq/PG_rep1.sorted.bam

# 2. GeneMark_ES Abinitio prediction
#perl ~/gm_et_linux_64/gmes_petap.pl --sequence ${ref} --ES --fungus --cores 14
#python3 ~/bin/gtfUtils.py -i genemark.gtf -reformat > genemark_ES.gtf
#awk 'BEGIN{FS=OFS="\t"}{gsub("\"","",$9);gsub("GeneMark.hmm","GeneMark_ES",$2);if($3~/gene/){gsub("gene_id ","ID=",$9)}else{if($3~/transcript/){gsub("transcript","mRNA", $3);gsub("gene_id ","Parent=",$9);gsub("transcript_id ","ID=",$9)}else{i=i+1;gsub("gene_id ","ID=",$9);gsub("_g","_g"i,$9)gsub("transcript_id ", "Parent=", $9)}};print}' genemark_ES.gtf > genemark_ES.gff3

# 3. GeneMark_ET integrate transcriptome
#~/gm_et_linux_64/star_to_gff.pl --star  SJ.out.tab --gff STAR.SJ2.gff --label STAR
#RNASEQ_hints=/disk/alpha/Morchella/gene_predict_STAR_hints/STAR.SJ2.gff
#/home/wanghm/gm_et_linux_64/gmes_petap.pl --sequence ${ref} --ET ${RNASEQ_hints} --fungus --cores 8
#python3 ~/bin/gtfUtils.py -i genemark.gtf -reformat > genemark_ES.gtf
#awk 'BEGIN{FS=OFS="\t"}{gsub("\"","",$9);gsub("GeneMark.hmm","GeneMark_ET",$2);if($3~/gene/){gsub("gene_id ","ID=",$9)}else{if($3~/transcript/){gsub("transcript","mRNA", $3);gsub("gene_id ","Parent=",$9);gsub("transcript_id ","ID=",$9)}else{i=i+1;gsub("gene_id ","ID=",$9);gsub("_g","_g"i,$9)gsub("transcript_id ", "Parent=", $9)}};print}' genemark_ET.gtf > genemark_ET.gff3

# 4. extract candicate coding regions from denovo transcripts assembly
# merge transcript assemblies
#stringtie --merge -o stringtie_assembly.merge.gtf -p 12 *.gtf
#trans_asm=/disk/alpha/OYR/stringtie_assembly.merge.gtf
#~/EvidenceModeler/TransDecoder-TransDecoder-v5.5.0/util/gtf_genome_to_cdna_fasta.pl ${trans_asm} ${ref} > transcripts.fasta
#/home/wanghm/EvidenceModeler/TransDecoder-TransDecoder-v5.5.0/util/gtf_to_alignment_gff3.pl ${trans_asm} > transcripts.gff3
#/home/wanghm/EvidenceModeler/TransDecoder-TransDecoder-v5.5.0/TransDecoder.LongOrfs -t transcripts.fasta
#/home/wanghm/EvidenceModeler/TransDecoder-TransDecoder-v5.5.0/TransDecoder.Predict -t transcripts.fasta
#/home/wanghm/EvidenceModeler/TransDecoder-TransDecoder-v5.5.0/util/cdna_alignment_orf_to_genome_orf.pl transcripts.fasta.transdecoder.gff3 transcripts.gff3 transcripts.fasta > transcripts.fasta.transdecoder.genome.gff3

########### run Evmodeler using transDecoder output; braker2 output;geneMark_ES output as gene_predictions.gff3
##########                using geneMark _ET output as transcript_alignments.gff3

#cat transcripts.fasta.transdecoder.genome.gff3 braker2.gff3 gene_predict_GeneMark_ES/genemark_ES.gff3 > gene_predictions.gff3
#ln gene_predict_STAR_hints/genemark_ET.gff3 transcript_alignments.gff3
#~/EvidenceModeler/EVidenceModeler-1.1.1/EvmUtils/partition_EVM_inputs.pl --genome ${ref} --gene_predictions gene_predictions.gff3 \
#	--transcript_alignments transcript_alignments.gff3 --segmentSize 100000 --overlapSize 10000 --partition_listing partitions_list.out
#~/EvidenceModeler/EVidenceModeler-1.1.1/EvmUtils/write_EVM_commands.pl --genome ${ref} --weights /disk/alpha/Morchella/evm_predict/weights.txt --gene_predictions gene_predictions.gff3 --transcript_alignments transcript_alignments.gff3 \
##				      --output_file_name evm.out  --partitions partitions_list.out >  commands.list
#~/EvidenceModeler/EVidenceModeler-1.1.1/EvmUtils/execute_EVM_commands.pl commands.list | tee run.log
#~/EvidenceModeler/EVidenceModeler-1.1.1/EvmUtils/recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out
#~/EvidenceModeler/EVidenceModeler-1.1.1/EvmUtils/convert_EVM_outputs_to_GFF3.pl  --partitions partitions_list.out --output evm.out  --genome ${ref}
#find . -regex ".*evm.out.gff3" -exec cat {} \; > EVM.all.gff3
# reformat GTF 
gffread -T EVM.all.gff3 -o EVM.all.gtf
python3 ~/bin/gtfUtils.py -i EVM.all.gtf -reformat -o EVM.all.reformat.gtf
python3 ~/bin/gtfUtils.py -i EVM.all.reformat.gtf -r S41_ -o S41.gtf
gffread -y S41.protein.fasta -T S41.gtf -g OYR.fasta
awk '/>/{split($2,a,"=");print ">"a[2];next}{print}' S41.protein.fasta > S41_protein.fasta
rm S41.protein.fasta
