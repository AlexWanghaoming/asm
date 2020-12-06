#!/bin/bash
set -eu
strain=CanS41

##################################### CAZymes
wget -c http://bcb.unl.edu/dbCAN2/data/blast/2020120505736/overview.txt -O CanS41_dbCan.txt
# get summary and Cazyme geneID
awk 'BEGIN{total["GH"]=0;total["GT"]=0;total["CE"]=0;total["AA"]=0;total["CBM"]=0;total["PL"]=0}FNR>1 && $6>=2{printf $1"\t";array["GH"]=0;array["GT"]=0;array["CE"]=0;array["AA"]=0;array["CBM"]=0;array["PL"]=0;for(i=2;i<=4;i++){for(k in array){if($i~k){array[k]++}}};for(k in array){if(array[k]>=2){printf k",";total[k]++}};printf "\n"}END{print "GH:"total["GH"]"\nGT:"total["GT"]"\nCE:"total["CE"]"\nAA:"total["AA"]"\nCBM:"total["CBM"]"\nPL:"total["PL"] > "CanS41_cazy.summary";close("CanS41_cazy.summary")}' CanS41_dbCan.txt > CanS41_cazy.geneID
awk 'NR==FNR{FS=" ";a[$1]=$1}NR>FNR{RS=">";FS="\n";split($0, id," ");if (id[1] in a) {printf ">%s",$0}}' CanS41_cazy.geneID CanS41_protein.fasta > CanS41_CAZYmes.fasta

################################## SM clusters
proxychains4 wget -c https://fungismash.secondarymetabolites.org/upload/fungi-b4e9e792-a67f-4953-a645-2a8ba083e37a/OYR.zip -O CanS41_antiSMASH.zip
unzip CanS41_antiSMASH.zip -d CanS41_antiSMASH
python3 antismash_gbk2fasta.py CanS41_antiSMASH/*.region*.gbk > CanS41_antiSMASHcoreGenes.fasta
awk -F ";" '/>/{a[$2]++}END{for(i in a){print i,a[i]}}' CanS41_antiSMASHcoreGenes.fasta > CanS41_antiSMASH.summary

########################### secretome
# deepLoc predictor
seqkit split -f -p 8 CanS41_protein.fasta # split fasta to 8 children
find CanS41_protein.fasta.split -name "*.fasta" | parallel -j 4 deeploc -f {} -o {}_deepLoc # 4 processing running
cat CanS41_protein.fasta.split/CanS41_protein.fasta.split/CanS41_protein*_deepLoc.txt > CanS41_deeploc.txt
awk '$2~/Extracellular/{print $1}' CanS41_deeploc.txt > CanS41_filterByDeeloc.geneID
rm -rf CanS41_protein.fasta.split/
# signalP predictor
signalp -fasta /home/wanghm/whm/pacbio_data/OYR/CanS41_protein.fasta -mature -gff3 -prefix ${strain}_signalp
# predict transmembrane helices
~/software/tmhmm-2.0c/bin/tmhmm ${strain}_signalp_mature.fasta > transme_helices.txt
# predict subcellular localization
~/software/targetp-2.0/bin/targetp -fasta ${strain}_signalp_mature.fasta -mature -gff3 -prefix ${strain}_targetP
# remove transmembrane helices
grep "Number of predicted TMHs:  0" transme_helices.txt | awk '{print $2}' | sort | uniq | grep -f - \
    ${strain}_signalp.gff3 > secretome_rmHelices.txt
# remove mitochondrial localization
awk '$0~/Mitochondrion/{print $1}' ${strain}_targetP.gff3 | grep -v -f - secretome_rmHelices.txt | awk -F "\t" '{print $1}' | sort | uniq > ${strain}_filterBySignalP.geneID
rm secretome_rmHelices.txt
rm transme_helices.txt
rm CanS41_signalp*
rm CanS41_targetP*
cat CanS41_filterBySignalP.geneID <(awk '{print $1}' CanS41_filterByDeeloc.geneID) | sort | uniq > CanS41_secretome.geneID
# extract secretome sequence
awk 'NR==FNR{a[$1]=$1}NR>FNR{RS=">";FS="\n";split($0, id," ");if (id[1] in a) {printf ">%s",$0}}' CanS41_secretome.geneID CanS41_protein.fasta > CanS41_secretome.fasta

################################ Effectors
python ~/software/EffectorP_2.0/Scripts/EffectorP.py -i CanS41_secretome.fasta -o CanS41_effector.txt -E CanS41_effector.fasta
grep "Effector probability" CanS41_effector.txt | awk -F "[ ]" '{print $1}' | sort | uniq > CanS41_effector.geneID

