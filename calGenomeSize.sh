awk 'BEGIN{i=1}/>/&&NR>1{print ""}{printf "%s",/^>/ ? i++" " : $0}' ../yn56.fasta | awk 'BEGIN{print "chr\tsize"}{print$1"\t"length($2)}' > genome/yn56.genome.txt
awk 'BEGIN{i=1}/>/&&NR>1{print ""}{printf "%s",/^>/ ? i++" " : $0}' ../hn47.fasta | awk 'BEGIN{print "chr\tsize"}{print$1"\t"length($2)}' > genome/hn47.genome.txt
awk 'BEGIN{i=1}/>/&&NR>1{print ""}{printf "%s",/^>/ ? i++" " : $0}' ../qz.fasta | awk 'BEGIN{print "chr\tsize"}{print$1"\t"length($2)}' > genome/qz.genome.txt
awk 'BEGIN{i=1}/>/&&NR>1{print ""}{printf "%s",/^>/ ? i++" " : $0}' ../fj11.fasta | awk 'BEGIN{print "chr\tsize"}{print$1"\t"length($2)}' > genome/fj11.genome.txt
awk 'BEGIN{i=1}/>/&&NR>1{print ""}{printf "%s",/^>/ ? i++" " : $0}' ../gz15.fasta | awk 'BEGIN{print "chr\tsize"}{print$1"\t"length($2)}' > genome/gz15.genome.txt
awk 'BEGIN{i=1}/>/&&NR>1{print ""}{printf "%s",/^>/ ? i++" " : $0}' ../gd10.fasta | awk 'BEGIN{print "chr\tsize"}{print$1"\t"length($2)}' > genome/gd10.genome.txt
