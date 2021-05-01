
#### YDJ circos ####
library(circlize)
library(Rsamtools)
library(GenomicFeatures)

## prepare data
ydj.fai <- read.table("/Users/alexwang/0data/Morchella/ydj.fasta.fai")
canu.gs =  ydj.fai$V2
txdb <- makeTxDbFromGFF("/Users/alexwang/0data/Morchella/SCLS.gtf", format = "gtf")
ge <- genes(txdb)

## GC content
GC = 0.4737
names(canu.gs) <- paste0("Scaffold_", seq(1,42))
ydj.window <- tileGenome(canu.gs, tilewidth = 10000, cut.last.tile.in.chrom = T)
bin_seq <- getSeq(FaFile("/Users/alexwang/0data/Morchella/ydj.fasta"), ydj.window)
gc_content <- letterFrequency(bin_seq, "GC") / letterFrequency(bin_seq, "ATCG")
ydj.window$gc <- gc_content-GC
## repeat 
repeat.df <- read.table("/Users/alexwang/0data/Morchella/ydj.fasta.out.gff")
## 6mA
AT_count <- letterFrequency(bin_seq, "AT")
m6AInWindos_count <- countOverlaps(ydj.window, m6AGrange)
ydj.window$m6A <- m6AInWindos_count / AT_count


chro <- ydj.fai$V1
starts = rep(0, length(canu.gs))
ends = canu.gs
genoCir <- data.frame(chr=chro, start=starts, end=ends)
genoCir$chr <- as.vector(genoCir[,1])
#### start plot######################################### only show fisrt 28 scaffolds of Morchella genome
circos.clear()
circos.par(start.degree = 87, track.height = 0.02, cell.padding = c(0,0,0,0), gap.degree=c(rep(1,26), 5))
circos.genomicInitialize(data = genoCir[1:27,],
                         sector.names = sub("Scaffold_","",chro[1:27]),
                         labels.cex = 0.5, track.height = 0.05, plotType = "labels")
### a: ideagram of 27 Scaffolds
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
    idx = which(chro[1:27] == sector.index)
    circos.axis(labels.cex = 0.5,direction = "outside", labels.niceFacing = T, labels = "", minor.ticks = 5, lwd = 0.8, 
                major.at = c(0, canu.gs[idx]), major.tick.length = 0.4)
  }, track.height = 0.02, bg.border = NA)

### b: plot gene heatmap
ge.df <- data.frame(chr = seqnames(ge), start= start(ge), end=end(ge))
## calculate gene density in slide windows
bed.gene <- genomicDensity(region = ge.df, window.size = 10000, overlap = F)
cc1 = colorRamp2(breaks = c(0, 1), colors = c("white","lightblue"))
circos.genomicTrackPlotRegion(bed.gene, track.height=0.05, bg.border=NA,
                              panel.fun = function(region, value, ...){
                                circos.genomicRect(region, value, col = cc1(value),border = NA, lwd=0.1)
                              })
### c: repeat heatmap
repeat.df <- data.frame(chr=repeat.df$V1, start=repeat.df$V4, end=repeat.df$V5)
bed.repeat <- genomicDensity(region = repeat.df, window.size = 10000, overlap = F)
cc2 <- colorRamp2(breaks = c(0,1), colors = c("white", "black"))
circos.genomicTrackPlotRegion(bed.repeat, track.height = 0.05, bg.border = NA,
                              panel.fun = function(region, value, ...){
                                  circos.genomicRect(region, value, col = cc2(value), border = NA, lwd = 0.1)
                              })
### d: CAZYmes distribution 
# cazy.geneID <- read.table("/Users/alexwang/0data/Morchella/cazy/ydj_cazy.geneID", header = T)[,1]
# cazy.gr <- ge[ge$gene_id %in% cazy.geneID]
# bed.cazy <- data.frame(chr=seqnames(cazy.gr), start=start(cazy.gr), end=end(cazy.gr))
# bed.cazy$value <- 1
# circos.genomicTrackPlotRegion(bed.cazy, track.height = 0.05,ylim = c(0,1), bg.border="gray",
#                               panel.fun = function(region, value, ...){
#                                 circos.genomicRect(region, value, col = "darkorange", border = "darkorange", lwd = 1)
#                               })

### d: 6mA heatmap
m6A.df <- data.frame(chr=seqnames(ydj.window), start=start(ydj.window), end=end(ydj.window), value = ydj.window$m6A)
# max(ydj.window$m6A)
# aaa <- boxplot(ydj.window$m6A)
cc3 <- colorRamp2(breaks = c(0,0.008), colors = c("white", "red"))
circos.genomicTrackPlotRegion(m6A.df, track.height = 0.05, bg.border = NA,
                              panel.fun = function(region, value, ...){
                                circos.genomicRect(region, value, col = cc3(value), border = NA, lwd = 0.1)
                              })
### e: full-length LTR 
ltr.df <- read.table( "/Users/alexwang/0data/Morchella/LTR/ydj.fasta.pass.list")
spl <- strsplit(ltr.df$V1, "[:.]")
chr = unlist(lapply(spl, function(x) x[1] ))
s = as.numeric(unlist(lapply(spl, function(x) x[2] )))
e = as.numeric(unlist(lapply(spl, function(x) x[4] )))
LTR.bed <- data.frame(chr = chr, start = s, end = e, type = ltr.df$V10)
LTR.bed$value <- 1
LTR.bed[LTR.bed$end < LTR.bed$start, ]
LTR.bed.list = list(LTR.bed[LTR.bed$type == "Gypsy", ], 
                    LTR.bed[LTR.bed$type == "Copia",  ],
                    LTR.bed[LTR.bed$type != "Gypsy" & LTR.bed$type != "Copia",  ])
circos.genomicTrackPlotRegion(LTR.bed.list, track.height = 0.05,ylim = c(0,1), bg.border="lightgray",
                              panel.fun = function(region, value, ...){
                              i=getI(...)
                              if(i == 1){
                                circos.genomicRect(region, value, col = "orange", border = "orange", lwd = 1)
                              } else if(i == 2){
                                circos.genomicRect(region, value, col = "red", border = "red", lwd = 1)
                              }else{
                                circos.genomicRect(region, value, col = "blue", border = "blue", lwd = 1)
                              } })


## f: SM clusters
sm.geneID <- read.table("/Users/alexwang/0data/Morchella/sm/scls/scls_SM_cluster.geneID")[,1]
sm.gr <- ge[ge$gene_id %in% sm.geneID]
bed.sm <- data.frame(chr=seqnames(sm.gr), start=start(sm.gr), end=end(sm.gr))
bed.sm$value <- 1
circos.genomicTrackPlotRegion(bed.sm, track.height = 0.05,ylim = c(0,1),bg.border = "lightgray",
                              panel.fun = function(region, value, ...){
                                circos.genomicRect(region, value, col = "darkgreen", border = "darkgreen", lwd = 1)
                              })

# custom_scale <- function(v, mi, ma){
#   (v - min(v)) * (ma-mi) / (max(v)-min(v)) + mi
# }
### g: RNAseq coverage barplot
RNA.df <- read.table("/Users/alexwang/0data/Morchella/ydj_RNAseq.regions.bed")
colnames(RNA.df) <- c("chr", "start", "end", "value")
RNA.df$value <- log2(RNA.df$value + 1)
circos.genomicTrackPlotRegion(RNA.df,track.height=0.08, bg.border=NA, ylim = c(min(RNA.df$value), max(RNA.df$value)),
                              panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, value, ybottom = 0, ytop.column = 1, 
                                                   col = "SlateBlue1", border = NA, lwd = 0.5, area = T)
                              })
### h: plot GC content
bed.gc <- data.frame(chr = seqnames(ydj.window), start = start(ydj.window), 
                     end = end(ydj.window))
vv <- ydj.window$gc
vv[which(vv %in% boxplot.stats(vv)$out)] <- mean(vv)
bed.gc$value <- vv
bed.gc.list <- list(bed.gc[bed.gc$value>0, ], 
                    bed.gc[bed.gc$value<0, ])
circos.genomicTrackPlotRegion(bed.gc.list, track.height = 0.08, bg.border = NA, ylim = c(min(vv),max(vv)),
                              # ylim=c(min(bed.gc$value), max(bed.gc$value)),
                              panel.fun = function(region, value, ...){
                                i=getI(...)
                                if(i == 1){
                                  # circos.genomicLines(region, value, col = "forestgreen", border = NA,
                                  #                     baseline = 0,area = T)
                                  circos.genomicRect(region, value, col = "forestgreen", border = NA,
                                                     ybottom = 0, ytop.column = 1)
                                }else{
                                  # circos.genomicLines(region, value, col = "firebrick1", border = NA,
                                  #                     baseline = 0,area = T)
                                  circos.genomicRect(region, value, col = "firebrick1", border = NA,
                                                     ybottom = 0, ytop.column = 1)
                                }})

##############################################
#### YDJ NZTD comparison circos ####
library(karyoploteR)
library(GenomicFeatures)
library(magrittr)
############ plot colin relationship between NZTD and SCLS
gs <- read.table(file = "/Users/alexwang/0data/Morchella/colin/nztd.genome.txt", header = T)
# gs <- read.table(file = "/Users/alexwang/0data/Morchella/ydj.genome.txt", header = T)
gs <- gs[order(gs$size),]
gs.gr <- GRanges(seqnames = gs$chr, ranges = IRanges(start=1, end = gs$size), col="red")
## bundle regions
max.gap <- 1000
bundle1 <- read.table("/Users/alexwang/0data/Morchella/colin/ydj.bundle_telo.txt", header = T)
bundle1 <- bundle1[abs(bundle1$s1-bundle1$e1)>max.gap & abs(bundle1$s2-bundle1$e2)>max.gap, ]

# sum(bundle1[bundle1$ref == 6, ]$e1 - bundle1[bundle1$ref == 6, ]$s1)

cc <- rainbow(60)[bundle1$ali]
x <- GRanges(seqnames = bundle1$ref, 
             ranges = IRanges(start = bundle1$s1, end = bundle1$e1),col = "red")
bundle.gr <- x
width(bundle.gr[seqnames(bundle.gr)=="4"])

pp = getDefaultPlotParams(plot.type = 1)
pp$leftmargin = 0.05
pp$topmargin = 80
pp$bottommargin = 15
pp$ideogramheight = 20
pp$data1inmargin = 10
pp$data1outmargin = 0
kp = plotKaryotype(genome = gs.gr, plot.params = pp, cex=0.6)
kpAddMainTitle(kp, "ydj_CoGe", cex=1)
kp <- kpPlotRegions(kp, data=bundle.gr,r0 = 0, r1 = 0.8,
                    border = NA, col = bundle.gr$col)

# telo1 <- bundle1$telo1
# telo2 <- bundle1$telo2
# for (idx in 1:NROW(bundle.gr)){
#   if (!is.na(telo1[idx])) {
#     kp <- kpPoints(karyoplot = kp, chr = seqnames(bundle.gr[idx]), x = as.numeric(telo1[idx]), 
#                    y=0.5, r0 = 0, r1 =0.8, ymax = 1, col="red", cex=0.4)
#   }
#   if (!is.na(telo2[idx])) {
#     kp <- kpPoints(karyoplot = kp, chr = seqnames(bundle.gr[idx]), x = as.numeric(telo2[idx]), 
#                    y=0.5, r0 =0 , r1 = 0.8, ymax = 1, col="red",cex = 0.4)
#   }
# }
###### calculate synetic prob
total.cov.core = 0
nucmer.coord <- read.table("/Users/alexwang/0data/pacbio/Morchella/colin/ydj.tab.txt", skip=4)
core.cov <- sum(nucmer.coord$V5)
total.cov.core = total.cov.core + core.cov

strain.genome <- read.table("/Users/alexwang/0data/pacbio/Morchella/colin/nztd.genome.txt", header = T)
core.sum <- sum(strain.genome$size)
total.cov.core/core.sum

############## circos plot 
library(circlize)
library(Rsamtools)
library(GenomicFeatures)
## prepare data
ydj.fai <- read.table("/Users/alexwang/0data/Morchella/ydj.fasta.fai")
NZTD.fai <- read.table("/Users/alexwang/0data/Morchella/colin/nztd.fa.fai")
FG.gs = NZTD.fai$V2
canu.gs =  ydj.fai$V2
txdb <- makeTxDbFromGFF("/Users/alexwang/0data/Morchella/ydj_primaryTrans.gtf", format = "gtf")
ge <- genes(txdb)

## GC content
GC = 0.4737
names(canu.gs) <- paste0("Scaffold_", seq(1,42))
ydj.window <- tileGenome(canu.gs, tilewidth = 10000, cut.last.tile.in.chrom = T)
bin_seq <- getSeq(FaFile("/Users/alexwang/0data/Morchella/ydj.fasta"), ydj.window)
gc_content <- letterFrequency(bin_seq, "GC") / letterFrequency(bin_seq, "ATCG")
ydj.window$gc <- gc_content-GC
## repeat 
repeat.df <- read.table("/Users/alexwang/0data/Morchella/ydj.fasta.out.gff")

canuVsRef.df <- read.table("/Users/alexwang/0data/Morchella/colin/ydj.tab.txt", skip = 4)
ref.iv <- canuVsRef.df[,c(14,1,2)]
# ref.iv$V14 <- sub("^", "chr", ref.iv$V14)
# ref.iv$V14 <- sub("chrMt", "Mt", ref.iv$V14)
canu.iv <- canuVsRef.df[,c(15,3,4)]
colnames(canu.iv) <- c("chr", "s", "e")
colnames(ref.iv) <- c("chr", "s", "e")
chro <- c(NZTD.fai$V1, ydj.fai$V1)

starts = rep(0, length(chro))
ends = c(FG.gs, canu.gs)
genoCir <- data.frame(chr=chro, start=starts, end=ends)
genoCir$chr <- as.vector(genoCir[,1])
## start plot
circos.clear()
circos.par(start.degree = 5, track.height = 1.0, cell.padding = c(0,0,0,0))
circos.genomicInitialize(data = genoCir, plotType = NULL)

sectors1 <- paste0("Scaffold_", seq(1,42))
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  if(sector.index %in% sectors1){
    idx=which(chro[60:101]==sector.index)
    chrom.tmp = chro[60:101]
    circos.text(mean(xlim), mean(ylim), labels =  sub("Scaffold_","",chrom.tmp[idx]), cex = 0.4, facing = "bending.inside", niceFacing = TRUE, family = "Sans" )
  }}, track.height = 0.02, bg.border = NA)

circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  if(sector.index %in% sectors1){
    idx = which(chro[60:101] == sector.index)
    circos.axis(labels.cex = 0.5,direction = "inside", labels.niceFacing = T, labels = "", minor.ticks = 0, lwd = 0.8, 
                major.at = c(0, canu.gs[idx]), major.tick.length = 0.4)
  }}, track.height = 0.02, bg.border = NA)
# chrom rect 
bed1 = data.frame(chr=paste0("Scaffold_", seq(1,42)), s = rep(1,42), e = canu.gs)
bed1$value <- sample(1:5, nrow(bed1), replace = T)
circos.genomicTrackPlotRegion(bed1, track.height = 0.02,bg.border=NA,ylim = c(0,1),
                              # track.margin=c(0.3,0.3),
                              panel.fun = function(region,value,...){
                                sector.index = get.cell.meta.data("sector.index")
                                if(sector.index %in% sectors1) {
                                  idx=which(sectors1==sector.index)
                                  circos.genomicRect(region,value,col = rainbow(42)[idx],border = NA)
                                }})
## plot gene density
bed.gene <- data.frame(chr = seqnames(ge), start = start(ge), end = end(ge))
circos.genomicDensity(bed.gene, track.height = 0.05, bg.border = NA, window.size = 10000)
# plot GC content
bed.gc <- data.frame(chr = seqnames(ydj.window), start = start(ydj.window), 
                     end = end(ydj.window))
bed.gc$value <- ydj.window$gc
circos.genomicTrackPlotRegion(bed.gc, track.height = 0.05, bg.border = NA, ylim=c(min(bed.gc$value), max(bed.gc$value)),
                              panel.fun = function(region, value, ...){
                                sector.index = get.cell.meta.data("sector.index")
                                if(sector.index %in% sectors1) {
                                  circos.genomicLines(region, value, col = "darkseagreen2", border = NA, area = F, lwd = 0.5)
                                }
                              })
# plot repeat regions
bed.repeat <- data.frame(chr=repeat.df$V1, start=repeat.df$V4, end=repeat.df$V5)
bed.repeat$value = 1
circos.genomicTrackPlotRegion(bed.repeat, track.height = 0.05, bg.border = NA, ylim=c(0,1),
                              panel.fun = function(region, value, ...){
                                sector.index = get.cell.meta.data("sector.index")
                                if(sector.index %in% sectors1) {
                                  circos.genomicRect(region, value, col = "darkorchid3", border = NA, lwd = 0.1)
                                }
                              })
########## inner track
sectors3 <- paste0("NZTD_", seq(1, 59))
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  if(sector.index %in% sectors3){
    idx=which(chro[1:59]==sector.index)
    chrom.tmp = chro[1:59]
    circos.text(mean(xlim), mean(ylim), labels =  sub("NZTD_", "" ,chrom.tmp[idx]), cex = 0.4, facing = "bending.inside", niceFacing = TRUE)
  }}, track.height = 0.02, bg.border = NA)

circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  if(sector.index %in% sectors3){
    idx = which(chro[1:59] == sector.index)
    circos.axis(labels.cex = 0.5,direction = "inside", labels.niceFacing = T, labels = "", minor.ticks = 0, lwd = 0.8, 
                major.at = c(0, FG.gs[idx]), major.tick.length = 0.4)
  }}, track.height = 0.02, bg.border = NA)

link.cc <- ifelse(ref.iv$chr == "NZTD_1", "red", "#BEBEBE80")
circos.genomicLink(ref.iv, canu.iv, border = NA, col =  link.cc, lwd = 0.5)
circos.clear()


###############################################
#### LTR distance result from LTR_retriever ####
library(ape)
library(Biostrings)
library(GenomicFeatures)
library(Rsamtools)
library(systemPipeR)

ltr <- "/Users/alexwang/0data/Morchella/LTR/ydj.fasta.pass.list"
genome <- "/Users/alexwang/0data/Morchella/ydj.fasta"

LTR_insert <- function(ltr, genome){
  ltr.df <- read.table(ltr)
  a1 <- strsplit(ltr.df$V1, split="[.:]", perl=T)
  chrom <- unlist(lapply(a1, function(x) x[1]))
  s_lLTR <- as.numeric(unlist(lapply(a1, function(x) x[2])))
  e_rLTR <- as.numeric(unlist(lapply(a1, function(x) x[4])))
  
  a2 <- strsplit(ltr.df$V7, split="[.:]", perl=T)
  e_lLTR <- as.numeric(unlist(lapply(a2, function(x) x[2]))) - 1
  s_rLTR <- as.numeric(unlist(lapply(a2, function(x) x[4]))) + 1
  
  genome.fasta <- Rsamtools::FaFile(genome)
  l.ltr <- GRanges(seqnames = chrom, 
                   ranges = IRanges(start = s_lLTR, end = e_lLTR))
  r.ltr <- GRanges(seqnames = chrom, 
                   ranges = IRanges(start = s_rLTR, end = e_rLTR))
  
  l.ltr.seq <- getSeq(genome.fasta, l.ltr)
  r.ltr.seq <- getSeq(genome.fasta, r.ltr)
  
  pair_align <- pairwiseAlignment(l.ltr.seq, r.ltr.seq)
  seq1 <- as.character(pattern(pair_align))
  seq2 <- as.character(subject(pair_align))
  res <- list()
  # res2 <- list()
  for (i in 1:length(seq1)) {
    x <- list(strsplit(seq1, "")[[i]], strsplit(seq2, "")[[i]])
    xx <- as.DNAbin(x)
    d <- dist.dna(xx)[1]
    # res2[[i]] <- d
    # fungal substitution rate  nucleotides per site per year
    r <- 1.05*10e-9 
    divergent_time <- d/(2*r)
    res[[i]] <- divergent_time/1000000
  }
  rr1 <- unlist(res)
  return(rr1)
}
ydj.ltr.res <- density(LTR_insert(ltr, genome))
options(scipen = 200)
p3 <- ~plot(ydj.ltr.res, main = "Time of insertion of LTR",col="purple", 
            lwd=1.5, xlab="Age (MYA)",ylim=c(0,4.2), xlim=c(0,3)
            # xlim=c(0,9), ylim=c(0, 0.3)
)

###### repeat pie plot 
retro = 2.21
transpoon = 9.15
unclassi = 4.72
simple = 0.94
pie.data <- data.frame(group = c("Retroelement", "DNA transposon", "Unclassified repeat", "Simple repeat"), 
                       ratio=c(retro, transpoon, unclassi, simple))
aa <- pie.data$ratio/sum(pie.data$ratio)
labs <- paste0(round(aa,3)*100,"%")
p <- ggpubr::ggpie(pie.data, x = "ratio", label = labs[c(1,2,3,4)], fill = "group", lab.pos = "in",
                   color="white", palette = RColorBrewer::brewer.pal(6, "Set2"), legend.title="")
pp <- ggpubr::ggpar(p, legend = "right", tickslab = F)
#### LTR class barplot 
gypsy <- 250
copia <- 35
unknow <- 5
par(mar=c(2, 5, 5, 2))
dd <- c(gypsy, copia, unknow)
names(dd) <- c("Gypsy", "Copia", "Unknow")
p2 <- ~barplot(dd, col = RColorBrewer::brewer.pal(3, "Set1"), ylab = "Counts")

###### LTR phylo tree
library(ggtree)
library(ggplot2)
library(treeio)
library(ape)
library(phytools)
tree <- read.tree("/Users/alexwang/0data/Morchella/LTR/LTR_internal.tree")
tree <- midpoint.root(tree)
LTR_class <- read.table("/Users/alexwang/0data/Morchella/LTR/LTR_class.txt")
colnames(LTR_class) <- c("label", "class")
LTR_class[LTR_class$class == "TCN1" | LTR_class$class == "MAGGY", ]$class = "Others"
##### sort tip.label by chrom and coord
# tiplabel.df <- stringr::str_split_fixed(tree$tip.label, pattern = "[_-]", n = 4)
# tiplabel.df <- as.data.frame(t(apply(tiplabel.df, 1, as.numeric)))
# order.idx <- with(tiplabel.df, order(V2,V3))
##### color tree by LTR insert time
# insert.age <- LTR_insert(ltr, genome)
# insertTime.df <- data.frame(label = tree$tip.label[order.idx],
#                             insert_time = insert.age)
# tree <- full_join(tree, insertTime.df, by="label")
# tree@phylo$tip.label <- paste0("LTR", 1:length(tree@phylo$tip.label))
# ggtree(tree, layout = "circular", lwd=0.2, aes(color=insert_time))  + 
#   scale_color_gradientn(colours=c("red","green",'blue'), limits=c(0,0.5)) + 
#   geom_tiplab2(size=0.6) 
#   # geom_text(aes(label=node), size=0.5)

# LTR_class$col <- rep(c("red","green"), 145)
# tree <- full_join(tree, LTR_class[,c(1,2)], by="label")
cls = split(LTR_class$label, LTR_class$class)
tree_OTU <- groupOTU(tree, cls)
p_tree <- ggtree(tree_OTU, aes(color=group), layout = "circular", branch.length = 'none') + 
  theme(legend.title = element_blank())

cowplot::plot_grid(plotlist = list(pp, p2, p3, p_tree), nrow = 2, ncol = 2, labels = "auto")

##############################################
#### plot FIR ####
# This script is used to plot intergenic regions length of secretome
## Generate Gene and Intergenic Ranges from GFF/GTF Annotations ## 
##  This is a modified version of the getFeat function by Thomas Girke
##	Computes intergenic, intron and overlapping gene ranges from GFF
##	or GTF files. Results can be combined with existing ranges. Input
##	and output objects are of class 'GRanges'.
getFeat2 <- function(x, gff, format="gff", range_types=c("intergenic", "gene_red", "intron", "gene", "exon")) {
  require(GenomicRanges)
  ## Input checks
  if(!missing(gff)) { x <- gff } # For backwards compatibility (if gff instead of x is used as function argument)
  if(class(x)!="GRanges") { stop("Input x needs of class GRanges (e.g. GFF/GTF file imported with import.gff)") }                                                                                                                                                                                                                   
  if(!format %in% c("gff", "GFF", "gtf", "GTF")) { stop("Format argument needs to be GFF or GTF") }
  
  ## If seqlengths (chr/space lengths) slot is empty, populate it with max value in each space (chr)
  if(any(is.na(seqlengths(x)))) {
    seqlengths(x) <- as.numeric(tapply(end(x), as.character(seqnames(x)), max))
  }
  
  ## Exon ranges
  elementMetadata(x) <- elementMetadata(x)[,c("type", "group", "score")]
  elementMetadata(x)[,"score"]<- strand(x)
  elementMetadata(x)[,"type"] <- as.character(as.data.frame(elementMetadata(x)[,"type"])[,1])
  if(format == "gff" | format == "GFF") { 
    elementMetadata(x)[,"group"] <- gsub("^.*?=(.*?)($|,|;|-).*", "\\1", elementMetadata(x)[,"group"])
  }
  if(format == "gtf" | format == "GTF") { 
    geneid <- gsub(";.*", "", elementMetadata(x)[,"group"])
    geneid <- gsub("\"| |gene_id", "", geneid)                                                                                                                                                       
    elementMetadata(x)[,"group"] <- geneid
  }
  strand(x) <- "*" # Erase strand information to make next steps strand insensitive
  x <- x[elementMetadata(x)[,"type"] != "chromosome", ] # Removes chromosome ranges from gff
  exon <- x[elementMetadata(x)[,"type"] == "exon",] # Returns only exon ranges
  
  ## Generate gene ranges when input is in GFF format
  if(format == "gff" | format == "GFF") { 
    gene <- x[elementMetadata(x)[,"type"] == "gene",] # Returns only gene ranges
    elementMetadata(gene)["length"]<- as.numeric(2 + end(gene) - start(gene))
  }
  
  ## Generate gene ranges when input is in GTF format
  if(format == "gtf" | format == "GTF") { 
    getgeneRanges <- function(x) {                                                                                                                                                                                                            
      tmpdf <- as.data.frame(x)                                                                                                                                                                                                              
      tmpgene <- tmpdf[order(tmpdf$group, tmpdf$start),]                                                                                                                                                                                       
      tmpgene <- tmpgene[!duplicated(tmpgene$group),]
      tmpgeneend <- tmpdf[order(tmpdf$group, -tmpdf$end),]                                                                                                                                                                                     
      tmpgeneend <- tmpgeneend[!duplicated(tmpgeneend$group),]                                                                                                                                                                                 
      tmpgene[,"end"] <- tmpgeneend$end                                                                                                                                                                                                        
      tmpgene[,"width"] <- tmpgene$end - tmpgene$start + 1                                                                                                                                                                                     
      tmpgene[,"type"] <- "gene"
      tmpgene[,"score"]<- tmpgene$score
      tmpgene <- tmpgene[order(tmpgene$seqnames, tmpgene$start),]                                                                                                                                                                              
      gr <- GRanges(seqnames = Rle(tmpgene$seqnames), ranges = IRanges(tmpgene$start, end = tmpgene$end), strand = Rle(tmpgene$strand))                                                                                                        
      elementMetadata(gr) <- tmpgene[,6:8]                                                                                                                                                                                                     
      seqlengths(gr) <- as.numeric(seqlengths(x))
      return(gr)                                                                                                                                                                                                                               
    }
    gene <- getgeneRanges(x)
    elementMetadata(gene)["length"]<- as.numeric(2 + end(gene) - start(gene))
  }
  
  ## Get range from chromosome start to first annotation feature                                         
  start <- tapply(start(gene), as.vector(seqnames(gene)), min)
  start <- paste(names(start), start, sep="_")
  all <- paste(as.vector(seqnames(gene)), start(gene), sep="_")
  firstinter <- gene[which(all %in% start),]
  firstinter <- unique(firstinter) 
  ranges(firstinter) <- IRanges(1, start(firstinter)-1)
  elementMetadata(firstinter)["group"] <- paste("start", as.character(as.data.frame(elementMetadata(firstinter)["group"])[,1]), sep="_")
  ## Get range from chromosome end to last annotation feature                 
  chrend <-  seqlengths(x) 
  end <- tapply(end(gene), as.vector(seqnames(gene)), max)
  end <- paste(names(end), end, sep="_")
  all <- paste(as.vector(seqnames(gene)), end(gene), sep="_")
  lastinter <- gene[which(all %in% end),]
  lastinter <- unique(lastinter) 
  ranges(lastinter) <- IRanges(end(lastinter)+1, chrend)
  elementMetadata(lastinter)["group"] <- paste(as.character(as.data.frame(elementMetadata(lastinter)["group"])[,1]), "end", sep="_")
  
  ## Reduced gene ranges
  if(any(range_types %in% c("gene_red", "intergenic"))) {
    ## Obtain gene ranges and collapse overlapping genes with reduce
    gene_red <- reduce(gene)
    ols <- as.matrix(findOverlaps(gene, gene_red))
    red_labels <- tapply(as.character(as.data.frame(elementMetadata(gene)["group"])[,1]), ols[,2], paste, collapse="_")
    elementMetadata(gene_red) <- data.frame(type="gene_red", group=red_labels)
    elementMetadata(gene_red)["type"] <- as.character(as.data.frame(elementMetadata(gene_red)["type"])[,1])
  }
  
  ## Intergenic ranges
  if(any(range_types %in% "intergenic")) {
    tmp_gene_red <- gene_red
    seqlengths(tmp_gene_red) <- rep(NA, length(seqlengths(tmp_gene_red)))
    intergenic <- gaps(tmp_gene_red)
    intergenic <- intergenic[start(intergenic) != 1]
    tmp <- intergenic; start(tmp) <- start(tmp)-1; end(tmp) <- end(tmp)+1
    ols <- as.matrix(findOverlaps(gene_red, tmp))
    inter_labels <- tapply(as.character(as.data.frame(elementMetadata(gene_red)["group"])[ols[,1],1]), ols[,2], paste, collapse="__")
    elementMetadata(intergenic) <- data.frame(type="intergenic", source=NA, phase=NA, group=inter_labels)
    metcol <- intersect(intersect(names(elementMetadata(firstinter)), names(elementMetadata(gene_red))), names(elementMetadata(lastinter)))
    intergenic <- c(firstinter[, metcol], intergenic[, metcol], lastinter[, metcol])
    intergenic <- intergenic[order(as.vector(seqnames(intergenic)), start(intergenic))]
    elementMetadata(intergenic)["type"] <- "intergenic"
    elementMetadata(intergenic)["length"]<- as.numeric(2 + end(intergenic) - start(intergenic))
  }
  
  ## Intron Ranges
  if(any(range_types %in% "intron")) {
    tmp_exon <- exon
    seqlengths(tmp_exon) <- rep(NA, length(seqlengths(tmp_exon)))
    exonRl <- split(tmp_exon, elementMetadata(tmp_exon)[,"group"])
    intron <- mendoapply(gaps, exonRl) 
    intron <- unlist(intron)
    intron <- intron[start(intron) != 1]
    elementMetadata(intron) <- data.frame(type="intron", group=names(intron))
    names(intron) <- NULL
    elementMetadata(intron)["type"] <- as.character(elementMetadata(intron)[,"type"])
    elementMetadata(intron)["group"] <- as.character(elementMetadata(intron)[,"group"])
  }
  
  ## Organinze data components in one GRanges object
  tmpall <- eval(parse(text=paste("c(", paste(range_types, collapse=", "), ")")))
  tmpall <- tmpall[order(as.vector(seqnames(tmpall)), start(tmpall))]
  return(tmpall)
}

filled.contour3 <-function (x = seq(0, 1, length.out = nrow(z)),
                            y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
                            ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
                            levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
                            col = color.palette(length(levels) - 1), plot.title, plot.axes, 
                            key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
                            axes = TRUE, frame.plot = axes,mar, ...) {
  # modification by Ian Taylor of the filled.contour function
  # to remove the key and facilitate overplotting with contour()
  # further modified by Carey McGilliard and Bridget Ferris
  # to allow multiple plots on one page
  
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    }
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")
  # mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  # on.exit(par(par.orig))
  # w <- (3 + mar.orig[2]) * par("csi") * 2.54
  # par(las = las)
  # mar <- mar.orig
  plot.new()
  # par(mar=mar)
  plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
  if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) 
    stop("no proper 'z' matrix specified")
  if (!is.double(z)) 
    storage.mode(z) <- "double"
  .filled.contour(as.double(x), as.double(y), z, as.double(levels), col = col)
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
  }
  else plot.axes
  if (frame.plot) 
    box()
  if (missing(plot.title)) 
    title(...)
  else plot.title
  invisible()
}

library(rtracklayer)
library(GenomicFeatures)

txdb <- makeTxDbFromGFF("/Users/alexwang/0data/Morchella/SCLS.gtf")
gene <- genes(txdb)
gff <- import.gff("/Users/alexwang/0data/Morchella/SCLS.gtf")
# gffgene <- getFeat2(x=gff, format="gtf", range_types=c("gene"))
# strand(gffgene)<-mcols(gffgene)$score
# mcols(gffgene)$score <- NULL
gffgene <- genes(txdb)
gffgene$length = width(gffgene)
gffgene$group <- gffgene$gene_id
# extract intergenic regions
gffintg <- getFeat2(x = gff,format = "gtf", range_types = "intergenic")
length_intg<- as.data.frame(cbind(seq(1:length(ranges(gffintg))), as.numeric(mcols(gffintg)$length)))
# mean(length_intg$V2)
colnames(length_intg)<-c("index", "length")
three_intg_index<-precede(gffgene, gffintg)
five_intg_index<-follow(gffgene, gffintg)
gene_data<- as.data.frame(cbind(as.character(mcols(gffgene)$group), as.character(strand(gffgene)), as.numeric(five_intg_index), as.numeric (three_intg_index)))
colnames(gene_data)<-c("geneid", "strand", "FivePrime_index", "ThreePrime_index")
tempdata<-merge(x=gene_data, y=length_intg, by.x="FivePrime_index", by.y="index",all.x=TRUE)
colnames(tempdata)<-c("delete1", "geneid", "strand", "ThreePrime_index", "fiveprime")

FIRdata<-merge(x=tempdata, y=length_intg, by.x="ThreePrime_index", by.y="index",all.x=TRUE)
FIRdata$ThreePrime_index<-NULL
FIRdata$delete1<-NULL
colnames(FIRdata)<- c("geneid", "strand", "fiveprime", "threeprime")

# bin break
NumBins = 40
if ((max(FIRdata$fiveprime, na.rm=TRUE)>max(FIRdata$threeprime, na.rm=TRUE)) == TRUE){
  FIR2Bin<-FIRdata$fiveprime} else {FIR2Bin<-FIRdata$threeprime}
FIR2Bin=FIR2Bin[which(FIR2Bin!=0)]
FIR2Bin<-na.omit(FIR2Bin)
BinSteps<-round(length(FIR2Bin)/(NumBins-1), digits=0)
FIR2BinOrd<-sort(FIR2Bin)
TempBinLimits<-FIR2BinOrd[seq(FIR2BinOrd[2*BinSteps], length(FIR2BinOrd),BinSteps)]
TempBinLimits[length(TempBinLimits)+1]<- max(FIR2Bin, na.rm=TRUE)
## 对窗口进行指数分布拟合
x<-seq(length(TempBinLimits))
fit<-nls(log(TempBinLimits) ~ a*x + b, start = c(a=0, b=0),algorithm = "port", weights = ((x-0.5*NumBins)^2))
pred=predict(fit, x)
BinLimits=c(1, round(exp(pred),0), max(FIR2Bin))

## Data binning
xbin=cut(FIRdata$fiveprime,breaks=c(BinLimits))
ybin=cut(FIRdata$threeprime, breaks= c(BinLimits))

FIRdata<-cbind(FIRdata, xbin, ybin, genevalue=rep(1,length(FIRdata$fiveprime)))
## calculate the exact position in each bin
kk <- function(f) {
  if(f %in% BinLimits){
    return(0)
  }else{
    aa = sort(c(BinLimits, f))
    lo.f <- (f - (aa[which(aa==f)+1]+aa[which(aa==f)-1])/2)/(aa[which(aa==f)+1]-aa[which(aa==f)-1])
    return(lo.f)
  }
}
x.offset <- unlist(apply(data.frame(FIRdata$fiveprime), 1, kk))
y.offset <- unlist(apply(data.frame(FIRdata$threeprime), 1, kk))
FIRdata$xoffset <- x.offset
FIRdata$yoffset <- y.offset

GenValMatrix<-with(FIRdata,tapply(genevalue, list(xbin, ybin), sum))

# plot Heatmap
x<-1:ncol(GenValMatrix)
y<-1:nrow(GenValMatrix)
zlim = range(as.numeric(unlist(GenValMatrix)) , finite=TRUE)

mypalette<-colorRampPalette(c( "white","darkblue", "forestgreen", "goldenrod1","orangered", "red3", "darkred"), space="rgb")
mycol=mypalette(2*max(GenValMatrix,na.rm=TRUE))
mylabels<-paste(BinLimits[1:length(BinLimits)-1],BinLimits[2:length(BinLimits)],sep="-", collapse=NULL)

filled.contour3(x, y, z=GenValMatrix, plot.title=title(main="", xlab = "Five prime intergenic regions",
                                                       ylab = "Three prime intergenic regions",
                                                       cex.main=1, cex.lab=1),
                key.title = title(main ="Number of genes", cex.main=0.5, line=1),col=mycol,
                levels=pretty(zlim,2*max(GenValMatrix,na.rm=TRUE)),
                plot.axes={
                  axis(1,at=x, labels=mylabels,las=2, cex.axis=0.7);
                  axis(2,at=y, labels=mylabels, las=2, cex.axis=0.7);
                }
)


##############################################
#### genome basic infomation plot ####
library(plotrix)
## subreads line density plot
# awk '!/>/{printf $0;next}{print ""; print $0}' subreads.fasta | awk '!/>/{print length($0)}' > subreads.length.txt
library(ggplot2)
subreads.len <- read.table("/Users/alexwang/0data/Morchella/subreads.length.txt", header = F, col.names = "length")
mapped.subreads.len <- read.table("/Users/alexwang/0data/Morchella/mapped_subreads.length", header = F, col.names = 'length')
reads.len <- data.frame(g = factor(c(rep(1, nrow(subreads.len)), rep(2, nrow(mapped.subreads.len)))), 
                        length = c(subreads.len$length, mapped.subreads.len$length))
p1 <- ggplot(reads.len, aes(x=length, fill=g)) + geom_density(alpha=0.5, outline.type = "full")+
  theme_classic() + xlab("Subreads length") + ylab("Subreads density (%)") + 
  scale_x_continuous(expand = c(0,0), limits = c(0,50000)) + scale_y_continuous(expand = c(0,0), labels = seq(0,0.0125, 0.0025)) + 
  scale_fill_discrete(guide=FALSE)
## exon numbers
library(GenomicFeatures)
txdb <- makeTxDbFromGFF("/Users/alexwang/0data/Morchella/SCLS.gtf", format = "gtf")
exon.gr.list <- exonsBy(txdb, by = "tx")
exon.num <- unlist(lapply(exon.gr.list, FUN = function(x) {NROW(x)}))

cc1 = c(colorRampPalette(c("#00b4d8", "#ade8f4"))(10), "#fb8500")  
exon.tbl <- table(exon.num)
# sum(exon.tbl[1:5])/length(exon.num)
exon.tbl[1]/length(exon.num)
exonMoreThan10 <- sum(exon.tbl[11:length(exon.tbl)])
exonNum.df <- data.frame(cc = c(as.vector(exon.tbl[1:10]), exonMoreThan10))
positions1 = rownames(exonNum.df)
p2 <- ggplot(data = exonNum.df, aes(x=rownames(exonNum.df),y=cc))+ geom_bar(stat = "identity",fill=cc1) + geom_text(label = exonNum.df$cc, nudge_y = 60) + 
  scale_x_discrete(limits=positions1, labels=c("1","2","3","4","5","6","7","8","9","10", ">10")) +
  scale_y_continuous(expand =c(0,0), limits = c(0,3600)) + theme_classic() + xlab("Exon number") + ylab("Number of genes")
# sum(exonNum.df$cc)
####### gene lengths
cc2 = c(colorRampPalette(c("#7b2cbf", "#e0aaff"))(8), "#538d22")
ge.length <- width(genes(txdb))
aa = table(cut(ge.length, breaks=seq(0,4000, 500)))
morethan4000 <- length(ge.length[ge.length > 4000])
geLen.df <- data.frame(cc = c(as.vector(aa), morethan4000))
x.lab = c("1-500", "501-1000", "1001-1500", "1501-2000", "2001-2500", "2501-3000", "3001-3500","3501-4000", ">4000")
positions2 = x.lab
p3 <- ggplot(data = geLen.df, aes(x = x.lab, y = cc)) + geom_bar(stat = "identity", fill=cc2) + geom_text(label = geLen.df$cc, nudge_y = 60) + 
  scale_x_discrete(limits = positions2, labels = positions2)  + 
  scale_y_continuous(expand = c(0,0),limits = c(0,3500)) + theme_classic() + xlab("Gene length distribution") + ylab("Number of genes")+
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0))
# sum(geLen.df$cc)

pp <- cowplot::plot_grid(plotlist = list(p1,p2,p3), align = "h", nrow = 1)
cowplot::save_plot("scls_genomeInfo.pdf", plot = pp, base_height = 5, base_width = 12)

## pacbio coverage
library(karyoploteR)
library(GenomicFeatures)
library(magrittr)
library(tidyverse)

gs <- read.table("/Users/alexwang/0data/Morchella/ydj.genome.txt", header = T)
gs.gr <- GRanges(seqnames = paste0("Contig_",gs$chr), ranges = IRanges(start=1, end = gs$size), col="red")
pp = getDefaultPlotParams(plot.type = 2)
pp$leftmargin = 0.05
pp$topmargin = 80
pp$bottommargin = 15
pp$ideogramheight = 20
pp$data1inmargin = 0
pp$data1outmargin = 0

kp = plotKaryotype(genome = gs.gr, plot.params = pp, cex=0.6)


# add track of minimum PacBio coverage in bins

## samtools depth ydj_pbSorted_mapped.bam > ydj_pb.depth
## python3 minBinDepth.py ydj_pb.depth ydj_pb.cov

bincov <- read.table("/Users/alexwang/0data/Morchella/ydj_pb.cov.txt")
bincov <- GRanges(seqnames = gsub("Scaffold","Contig",as.character(bincov$V1)),
                  ranges = IRanges(start=bincov$V2, end = bincov$V3),
                  value = as.numeric(bincov$V6))

bincov$value[bincov$value>100] = 100
# kp = kpLines(kp, chr = seqnames(bincov),
#              x = start(bincov) + (end(bincov) - start(bincov)) / 2, 
#              y = bincov$value, 
#              ymin = 0, ymax = 200,
#              col = "darkblue", lwd = 1.5, clipping = T, r0 = 0, r1 = 1)

kp <- kpArea(kp, chr = seqnames(bincov), x = start(bincov) + (end(bincov) - start(bincov)) / 2,
             y = bincov$value, ymin = 0, ymax = 100, col = "lightblue", border = F)

telo <- read.table("/Users/alexwang/0data/Morchella/ydj.telo.txt", header = T)[,c(1,3,4)]
# left
kp = kpPoints(kp, chr = paste0("Contig_", telo$ref), 
               x = as.numeric(telo$left), 
               y = 0, col = "red")
# right
kp = kpPoints(kp, chr = paste0("Contig_", telo$ref), 
              x = as.numeric(telo$right), 
              y = 0, col = "red")

##################################################
#### 6mA ####
library(GenomicFeatures)
library(Rsamtools)
library(RColorBrewer)
library(rtracklayer)
txdb <-makeTxDbFromGFF(file = "/Users/alexwang/0data/Morchella/SCLS.gtf") 
m6A_gff <-data.table::fread("/Users/alexwang/0data/Morchella/6mA/ydj_m6A.gff")
m6AGrange <- GRanges(seqnames = m6A_gff$V1,ranges = IRanges(start = m6A_gff$V4, end = m6A_gff$V5),strand = m6A_gff$V7)
m6AGrange_total <- GRanges(seqnames = m6A_gff$V1,ranges = IRanges(start = m6A_gff$V4, end = m6A_gff$V5),strand = "*") # m6A_total do not distinct pos/neg strands

library(seqLogo)
library(ggseqlogo)
m6A_context <- readDNAStringSet(filepath = "/Users/alexwang/0data/Morchella/6mA/m6A_context.fasta", format = "fasta")
matrixForSeqLogo1 <- consensusMatrix(DNAStringSet(m6A_context, start = 2, end = 9), as.prob = T)[1:4,]
# p1 <- seqLogo(matrixForSeqLogo1, ic.scale = T)
p1 <- ggseqlogo(matrixForSeqLogo1)  ## ggplot obj

comparedSeqLogo <- function(methylationFile){
  library(seqLogo)
  library(ggseqlogo)
  library(Biostrings)
  ph1_m6A <- readDNAStringSet(filepath = methylationFile, format = "fasta")
  matrixForSeqLogo1 <- consensusMatrix(DNAStringSet(ph1_m6A,start = 2, end = 9),as.prob = T)[1:4,]   # cut -2 +4 motifs
  ph1_genome <- readDNAStringSet(filepath = "/Users/alexwang/0data/Morchella/ydj.fasta",format = "fasta")
  merged_genome <- Reduce(c,ph1_genome)
  samp <- sample(matchPattern("NNNNNANNNNN",fixed = F,subject = merged_genome),length(ph1_m6A)) # fixed = False can match N
  matrixForSeqLogo2 <- consensusMatrix(samp, as.prob = T)[1:4,]
  p1 <- ggseqlogo(matrixForSeqLogo1)
  p2 <- ggseqlogo(matrixForSeqLogo2)
  gridExtra::grid.arrange(p1, p2)
  # library(seqLogo)
  # ss <- seqLogo(matrixForSeqLogo2,ic.scale = T)
}
comparedSeqLogo(methylationFile = "/Users/alexwang/0data/Morchella/6mA/m6A_context.fasta")

library(tidyverse)
library(GenomicFeatures)
library(Rsamtools)
library(RColorBrewer)
library(rtracklayer)
library(ggplot2)

### ************** count 6mA in Genes
ge <- genes(txdb)
siteCountsOngenes <- countOverlaps(ge,m6AGrange_total)
table(siteCountsOngenes)
# concentrate on methylated genes
target_gene.list <- names(siteCountsOngenes[siteCountsOngenes!=0])  # 15598 genes totally 
## plot barplot for 6mA-marked genes counts, 
p2 <- ggpubr::ggbarplot(reshape2::melt(data.frame(1425, NROW(siteCountsOngenes)-1425)), x="variable", y= "value", label = T, fill = "lightblue", width = 0.6)+
  scale_y_continuous(expand = c(0,0),limits = c(0,17000)) +
  scale_x_discrete(breaks = c("X1425", "NROW.siteCountsOngenes....1425"), labels = c("non-6mA genes","6mA_marked genes")) + 
  xlab("") + ylab("Number of genes") + theme(axis.text.x = element_text(vjust = -1,size = 10))

target.geneID <- names(siteCountsOngenes[siteCountsOngenes>10])
write.table(target.geneID, "/Users/alexwang/0data/Morchella/m6A_target.geneID", quote = F, row.names = F, col.names = F)

########### GO enricherment analysis
library(GO.db)
library(tidyverse)
# write GO base infomation
godb <- AnnotationDbi::select(GO.db, keys(GO.db), columns(GO.db))
# write_tsv(godb, "/Users/alexwang/0data/db/godb.txt")
interpro <- read.table("/Users/alexwang/0data/Morchella/scls_interpro.tsv", fill = T, sep = "\t", 
                       col.names = paste0("V", 1:15), na.strings = "")
# write gene2GO 
gene2go <- interpro %>% dplyr::select(Gene = V1, GOID = V14) %>%
  mutate(Gene = str_replace(Gene, "\\.\\d+$", "")) %>% na.omit() %>%
  separate(GOID, paste0("X", 1:(max(str_count(.$GOID,"\\|"))+1), seq = ""), sep = "\\|") %>%
  gather(key = "X", value = "GOID", -Gene) %>% dplyr::select(Gene, GOID) %>%
  na.omit() %>% base::unique()
gene2go$Gene <- gsub(".t.*", "",gene2go$Gene)

write_tsv(gene2go, "/Users/alexwang/0data/Morchella/gene2go.txt")
go_anno <- gene2go %>% left_join(godb)
write_tsv(go_anno, paste0("/Users/alexwang/0data/Morchella/go_annot.txt"))

library(clusterProfiler)
godb <- read.delim("/Users/alexwang/0data/db/godb.txt",header = T,col.names = c("term","def", "ontology", "name"))
term2name <- godb[,c(1,4)]
term2gene <- read.table("/Users/alexwang/0data/Morchella/gene2GoID.txt", header = T, col.names = c("gene", "term"), fill = T)[,c(2,1)]
er <- enricher(target.geneID, TERM2GENE = term2gene, TERM2NAME = term2name)
enricher()
clusterProfiler::dotplot(er)
clusterProfiler::cnetplot(er)

## RNAseq expression of genes with >10 6mA sits
tpm.df <- read.table("/Users/alexwang/0data/Morchella/scls_tpm.tsv")
x1 <- tpm.df[rownames(tpm.df) %in% target.geneID, ][,2]
x2 <- tpm.df[!rownames(tpm.df) %in% target.geneID, ][,2]
t.test(x1, x2, alternative = "greater")
boxplot(x1, x2, ylim = c(0, 50))
## ********************************************************************* ##
library(ChIPseeker)
library(tidyverse)
library(tidyselect)
library(ggpubr)
library(annotatePeak)

m6A_anno.gr <- annotatePeak(m6AGrange_total, tssRegion = c(-1000,-300), TxDb = txdb, 
                            genomicAnnotationPriority=c("Promoter","Exon", "Intron", "Intergenic"), 
                            addFlankGeneInfo = T) %>% as.GRanges()
pm <- promoters(txdb, upstream = 500, downstream = 300)
genes(txdb)
m6A_anno <- m6A_anno.gr$annotation

# vars selection
a1 <- length(starts_with("Exon", vars = m6A_anno))
a2 <- length(starts_with("Intron", vars = m6A_anno))
a3 <- length(starts_with("Promoter", vars = m6A_anno))
a6 <- length(starts_with("Downstream", vars = m6A_anno))
a7 <- length(starts_with("Distal Intergenic",vars = m6A_anno))

## pie plot
pie.data <- data.frame(group = c("Promoter", "Exon", "Intron", "Intergenic"), 
                       counts=c(a3,a1,a2,a6+a7))
aa <- pie.data$counts/sum(pie.data$counts)
labs <- paste0(round(aa,3)*100,"%")
p <- ggpubr::ggpie(pie.data, x = "counts", label = labs[c(1,2,3,4)], fill = "group", lab.pos = "in",
                   color="white", palette = RColorBrewer::brewer.pal(6, "Set2"), legend.title="")
p3 <- ggpar(p, legend = "right", tickslab = F)

## *******************************  
gs <- "/Users/alexwang/0data/Morchella/ydj.genome.txt"
gs <- read.table(gs,header = T)
gs.gr <- GRanges(seqnames = paste0("Scaffold_", gs$chr), IRanges(start = 1, end = gs$size))

calculate_feature_methyFra <- function(x) {
  allcdsSeq <- getSeq(FaFile("/Users/alexwang/0data/Morchella/ydj.fasta"), x)
  methylationFraction <- sum(countOverlaps(x, m6AGrange_total)) / sum(letterFrequency(allcdsSeq, letters = c("A", "T")))
  return(methylationFraction)
}  # calculate methylation ratio on gr obj

get_feature_methylation_ratio <- function(feature){
  # feature <- feature[names(feature) %in% target_gene.list]
  allgeneSeq <- getSeq(FaFile("/Users/alexwang/0data/Morchella/ydj.fasta"), feature)
  ph1_gene <- feature # feature is GRanges Object
  
  # count each gene binding methylation sites
  allgeneCounts <- countOverlaps(ph1_gene, m6AGrange_total)
  
  # get each gene AT counts
  A_counts_perGenes <- letterFrequency(allgeneSeq, letters = "A")
  T_counts_perGenes <- letterFrequency(allgeneSeq, letters = "T")
  # total methylation ratio
  total_m6AFra <- allgeneCounts / (A_counts_perGenes + T_counts_perGenes)
  
  # calculate strand-specific methylation ratio
  sameDir_withGenes_6mACounts <- countOverlaps(ph1_gene, m6AGrange)
  diffDir_withGenes_6mACounts <- allgeneCounts - sameDir_withGenes_6mACounts
  samDir_6mAFra <- sameDir_withGenes_6mACounts / A_counts_perGenes
  diffDir_6mAFra <- diffDir_withGenes_6mACounts / T_counts_perGenes
  
  # make df merged
  m6A_fra_df <- data.frame(total=total_m6AFra, sameDir=samDir_6mAFra, diffDir=diffDir_6mAFra, row.names = names(ph1_gene))
  colnames(m6A_fra_df) <- c("total","sameDir","diffDir")
  # m6A_fra_df <- m6A_fra_df[m6A_fra_df[,1]!=0,]
  # m6A_fra_df <- na.omit(m6A_fra_df)
  return(m6A_fra_df)
}   # *** strand_specific methylation ratio

pmt.gr <- promoters(genes(txdb), upstream = 1000)
pmt.gr <- resize(pmt.gr, 700, fix = "start", ignore.strand = F)
pmt.gr <- subsetByOverlaps(pmt.gr, gs.gr, type = "within")
cds.gr <- cds(txdb)
exon.gr <- exons(txdb)
intron.gr <- unlist(intronsByTranscript(txdb), use.names = F)
## extract genic and intergenic region as GRanges object 
genic.gr <- GenomicRanges::reduce(ge, ignore.strand=T)
genic.gr <- genic.gr[start(genic.gr) > 0]
genic.gr <- subsetByOverlaps(genic.gr, gs.gr, type = "within")
intergenic <- gaps(genic.gr)
intergenic.gr <- intergenic[strand(intergenic)=="*"]

b1 <- calculate_feature_methyFra(pmt.gr)
b2 <- calculate_feature_methyFra(cds.gr)
b3 <- calculate_feature_methyFra(exon.gr)
b6 <- calculate_feature_methyFra(intron.gr)
b8 <- calculate_feature_methyFra(intergenic.gr)

methyRatio <- data.frame(type=as.factor(c("Promoter","Exon", "Intron", "Intergenic")), 
                         ratio=c(b1,b3,b6,b8)*100)
## plot feature's methylation ratio
 ggplot() + geom_bar(data = methyRatio, aes(x=type, y=ratio, fill=type), stat = "identity") + 
  scale_y_continuous(expand =c(0,0), limits = c(0,0.6)) + xlab("") + ylab("6mA / AT (%) ") + 
  theme_classic() + guides(fill=F) + scale_fill_manual(values = RColorBrewer::brewer.pal(4, "Set2"))

cowplot::plot_grid(plotlist = list(p1,p2,p3,p4), nrow = 2, align = c("h","v"), labels = "auto")

############################# gene 6mA ratio
library(reshape2)
get_feature_methylation_ratio <- function(feature){
  allgeneSeq <- getSeq(FaFile("/Users/alexwang/0data/Morchella/ydj.fasta"), feature)
  ph1_gene <- feature # feature is GRanges Object
  
  # count each gene binding methylation sites
  allgeneCounts <- countOverlaps(ph1_gene, m6AGrange_total)
  
  # get each gene AT counts
  A_counts_perGenes <- letterFrequency(allgeneSeq, letters = "A")
  T_counts_perGenes <- letterFrequency(allgeneSeq, letters = "T")
  # total methylation ratio
  total_m6AFra <- allgeneCounts / (A_counts_perGenes + T_counts_perGenes)
  
  # calculate strand-specific methylation ratio
  sameDir_withGenes_6mACounts <- countOverlaps(ph1_gene, m6AGrange)
  diffDir_withGenes_6mACounts <- allgeneCounts - sameDir_withGenes_6mACounts
  samDir_6mAFra <- sameDir_withGenes_6mACounts / A_counts_perGenes
  diffDir_6mAFra <- diffDir_withGenes_6mACounts / T_counts_perGenes
  
  # make df merged
  m6A_fra_df <- data.frame(total=total_m6AFra, sameDir=samDir_6mAFra, diffDir=diffDir_6mAFra, row.names = names(ph1_gene))
  colnames(m6A_fra_df) <- c("total","sameDir","diffDir")
  return(m6A_fra_df)
}
m6A_fra_df <- get_feature_methylation_ratio(feature = genes(txdb))
le <- read.table("/Users/alexwang/0data/Morchella/ydj_exp.txt", header = T)
ydj.exp <- dplyr::mutate(ydj.exp, TPM = (ydj1+ydj2)/2)
merge.df <- merge(ydj.exp, m6A_fra_df, by.x="geneID", by.y= "row.names")
merge.df <- merge.df[merge.df$TPM > 0 & merge.df$total > 0, ]
merge.df <- merge.df[order(merge.df$TPM, decreasing = T),]
### split to 10 groups
tpm_level.list <- split(merge.df, (seq(nrow(merge.df))-1) %/% (nrow(merge.df)/10))
lst10 <- lapply(tpm_level.list, function(x) x[,5])
reshapped <- melt(lst10)
ggplot(data = reshapped, aes(x=L1, y=value, fill=L1)) + geom_boxplot() + ylim(c(0, 0.02))
### split to 3 groups
# high <- merge.df[merge.df$TPM>30, ]$total
# media <- merge.df[merge.df$TPM < 30 & merge.df$TPM>1, ]$total
# low <- merge.df[merge.df$TPM<1, ]$total
# dd <- data.frame(high=high, media=media, low=low)
# boxplot(list(high, media, low))

## 6mA with repeat
library(regioneR)
library(ggsignif)
ltr <- "/Users/alexwang/0data/Morchella/LTR/ydj.fasta.pass.list"
ltr.df <- read.table(ltr)
spt1 <- strsplit(ltr.df$V1, split="[:.]")
chrom <- unlist(lapply(spt1, function(x) x[1]))
full_LTR.gr <- GRanges(seqnames = chrom, ranges = IRanges(start = as.numeric(unlist(lapply(spt1, function(x) x[2]))),
                                                          end = as.numeric(unlist(lapply(spt1, function(x) x[4])))))
spt2 <- strsplit(ltr.df$V7, split="[.:]", perl=T)
e_lLTR <- as.numeric(unlist(lapply(spt2, function(x) x[2]))) - 1
s_rLTR <- as.numeric(unlist(lapply(spt2, function(x) x[4]))) + 1
internal.gr <- GRanges(seqnames = chrom, ranges = IRanges(start=e_lLTR+1, end=s_rLTR-1))
t1 <- get_feature_methylation_ratio(internal.gr)[,1]
# dd = data.frame(type = factor(rep(c("LTR","Random regions"), each=nrow(ltr.df))),
#                 ratio = c(t1,t2))
# ggplot(dd, aes(x=type, y = ratio, fill=type)) + geom_boxplot() + 
#   theme_classic() + scale_y_continuous(expand = c(0,0), limits = c(0,0.06)) + xlab("") + ylab("6mA methylation ratio") + 
#   scale_fill_discrete(guide=FALSE) + 
#   geom_signif(aes(y_position=0.057),comparisons = list(c("LTR", "Random regions")),
#               test = "wilcox.test", test.args = "greater",tip_length=0.01)

repeat.df <- read.table("/Users/alexwang/0data/Morchella/ydj.fasta.out.gff")
repeat.gr <- GRanges(seqnames = repeat.df$V1, ranges = IRanges(start = repeat.df$V4, end = repeat.df$V5))
other_repeat.gr <- subsetByOverlaps(repeat.gr, full_LTR.gr, invert = T)
t2 <- get_feature_methylation_ratio(other_repeat.gr)[,1]

random.gr <- regioneR::createRandomRegions(nregions = nrow(ltr.df), genome = gs.gr)
t3 <- get_feature_methylation_ratio(random.gr)[,1]
random.gr <- regioneR::createRandomRegions(nregions = NROW(other_repeat.gr), genome = gs.gr)
t3 <- get_feature_methylation_ratio(random.gr)[,1]

marksig <- function(x1, x2, y, p, seg_len, offset) {
  lines(c(x1, x2), c(y, y), lwd = 2)
  lines(c(x1, x1), c(y, y - seg_len), lwd = 2)
  lines(c(x2, x2), c(y, y - seg_len), lwd = 2)
  text(x1+(x2-x1)/2, y + offset, p, cex = 1)
}
t.test(t1,t3)
par(mar=c(5,5,3,3), xpd = T)
boxplot(list(LTR=t1, Repeat=t2, Random=t3), ylim = c(0,0.025), 
        col = c("#00BFFF", "#FF7F50", "lightgoldenrod"), 
        axes = T,
        outline = F,
        frame = F,
        xlab = "", 
        ylab = "6mA methylation ratio (%)")

marksig(1,3, 0.025, "p-value < 2.2e-16", seg_len = 0.001, offset = 0.001)
marksig(2,3, 0.022, "p-value = 0.31", seg_len = 0.001, offset = 0.001)

#### 6mA profile ####
genome <- "/Users/alexwang/0data/Morchella/ydj.fasta"
m6A_gff <- fread("/Users/alexwang/0data/Morchella/6mA/ydj_m6A.gff", skip = 1)
m4C_gff <- fread("/Users/alexwang/0data/Morchella/6mA/ydj_m4C.gff", skip = 1)

m6AGrange = GRanges(seqnames = m6A_gff$V1,ranges = IRanges(start = m6A_gff$V4, end = m6A_gff$V5),strand = m6A_gff$V7)
m4CGrange = GRanges(seqnames = m4C_gff$V1,ranges = IRanges(start = m4C_gff$V4, end = m4C_gff$V5),strand = m4C_gff$V7)
m6AGrange_igstr <- m6AGrange
m4CGrange_igstr <- m4CGrange
strand(m6AGrange_igstr) = "*"
strand(m4CGrange_igstr) = "*"
m6AGrange2 <- m6AGrange
strand(m6AGrange2[strand(m6AGrange2)=="+"]) = "*"
strand(m6AGrange2[strand(m6AGrange2)=="-"]) = "+"
strand(m6AGrange2[strand(m6AGrange2)=="*"]) = "-"
###### get promoters
gs <- "/Users/alexwang/0data/Morchella/ydj.genome.txt"
gs <- read.table(gs,header = T)
gs.gr <- GRanges(seqnames = paste0("Scaffold_", gs$chr), IRanges(start = 1, end = gs$size))
txdb <-makeTxDbFromGFF(file = "/Users/alexwang/0data/Morchella/EVM.all.reformat.rename.gtf") 
promoter <- promoters(genes(txdb), upstream = 500, downstream = 1000)
## limit the promoter to genome size
promoter <- subsetByOverlaps(promoter, gs.gr, type = "within")

forEachStrand <- function(y) {### operate Views list
  tagMatrixList <-lapply(y, function(x) t(viewApply(x, as.vector)))
  tagMatrix <- do.call("rbind", tagMatrixList)
  rownames(tagMatrix) <- 1:nrow(tagMatrix)
  tagMatrix <- tagMatrix[rowSums(tagMatrix) != 0,]   # remove no methylation genes
  return(tagMatrix)
}

par(mfcol=c(2,1),xpd=T)
# par()
#### do not distinguish coding and template strand
getAroundTssMatrix_profile <- function(window,GRangeObj){
  peak.cov1 <- coverage(GRangeObj) 
  cov.len1 <- elementNROWS(peak.cov1)
  cov.width1 <- GRanges(seqnames = names(cov.len1), IRanges(start = rep(1,length(cov.len1)), end = cov.len1))
  windows1 <- subsetByOverlaps(window, cov.width1, type = "within")  # reconfirm all promoter regions are in genome
  chr.idx1 <- intersect(names(peak.cov1), unique(as.character(seqnames(windows1))))
  peakView1 <- Views(peak.cov1[chr.idx1], as(windows1[strand(windows1)=="+"], "IntegerRangesList")[chr.idx1])## showMethods("coerce") can get all avaliable transforms
  peakView2 <- Views(peak.cov1[chr.idx1], as(windows1[strand(windows1)=="-"], "IntegerRangesList")[chr.idx1])## showMethods("coerce") can get all avaliable transforms
  
  peakView1 <- peakView1[lapply(peakView1, length)>0]
  peakView2 <- peakView2[lapply(peakView2, length)>0]
  
  diff_strandList <- list(peakView1,peakView2)
  diffStrandMatrixList <- lapply(diff_strandList, forEachStrand)
  return(diffStrandMatrixList)
}
aroundTssProfile <- function(genomeFile, upstream = NULL,downstream = NULL, windowSize) {
  k=1
  for (i in 1:2) {
    ifelse(k==1,geneMethylationCount <- countOverlaps(promoter,m6AGrange), geneMethylationCount <- countOverlaps(promoter,m4CGrange))
    promoter_trimed <- promoter[names(promoter) %in% names(geneMethylationCount[geneMethylationCount!=0])]
    promoterSeq <- getSeq(FaFile(genomeFile), promoter_trimed) 
    ifelse(k==1,tt <- getAroundTssMatrix_profile(promoter, m6AGrange_igstr), tt <- getAroundTssMatrix_profile(promoter, m4CGrange_igstr))
    tt[[2]] <- t(apply(tt[[2]],1,rev))
    tt1 <- do.call(rbind,tt)
    ss <- colSums(tt1)
    ifelse(k==1,countBase <- (consensusMatrix(promoterSeq)[1,]+consensusMatrix(promoterSeq)[4,]), 
           countBase <- (consensusMatrix(promoterSeq)[2,]+consensusMatrix(promoterSeq)[3,])) 
    # promoterSeq <- getSeq(FaFile("/home/wanghm/whm/1FG/35genome/Fusarium_graminearum.RR1_with_mito.fa"),promoter_trimed)
    ss <- ss / countBase
    ss[which(is.nan(ss))] <- 0    # because YDJ gene annotations have no UTR, so 503 positions are all ATG
    pos <- value <- NULL
    tagCount <- data.frame(pos = c(-499:1000), value = ss)
    winPos <- matrix(tagCount$pos, ncol = windowSize, byrow = TRUE)
    winValue <- matrix(tagCount$value, ncol = windowSize, byrow = TRUE)
    xx = rowMeans(winPos)
    yy = rowMeans(winValue)*100
    if(k==1){    # ,ylim=c(0.0030, 0.0060)
      plot(xx,yy, type="l",col="#157a6e", xlab = "Distance from TSS(bp)",ylim=c(0.40, 0.60),ylab = paste("methylation occupancy (%)"), lwd=2)
      legend("topright", c("6mA"), lty = c(1,1,2), lwd = 2, 
             col = c("#157a6e"), bty = "n", cex = 0.8, xpd = TRUE)
      
    }
    # else{
    #   lines(xx,yy,lwd=2,col="#e2711d")
    #   
    # }
    k=k+1
  }
  
}
aroundTssProfile(genomeFile = "/Users/alexwang/0data/Morchella/ydj.fasta",windowSize = 50)

##### strand specific TSS profile
getAroundTssMatrix <- function(window,GRangeObj){
  peak.cov1 <- coverage(GRangeObj[strand(GRangeObj) == "+"]) # coverage return Rle Object
  peak.cov2 <- coverage(GRangeObj[strand(GRangeObj) == "-"]) # coverage return Rle Object
  cov.len1 <- elementNROWS(peak.cov1)
  cov.width1 <- GRanges(seqnames = names(cov.len1), IRanges(start = rep(1,length(cov.len1)), end = cov.len1))
  windows1 <- subsetByOverlaps(window, cov.width1, type = "within")  # reconfirm all promoter regions are in genome
  chr.idx1 <- intersect(names(peak.cov1), unique(as.character(seqnames(windows1))))
  peakView1 <- Views(peak.cov1[chr.idx1], as(windows1[strand(windows1)=="+"], "IntegerRangesList")[chr.idx1])## showMethods("coerce") can get all avaliable transforms
  
  
  # For peak.cov2
  cov.len2 <- elementNROWS(peak.cov2)
  cov.width2 <- GRanges(seqnames = names(cov.len2), IRanges(start = rep(1,length(cov.len2)), end = cov.len2))
  windows2 <- subsetByOverlaps(window, cov.width2, type = "within")  # reconfirm all promoter regions are in genome
  chr.idx2 <- intersect(names(peak.cov2), unique(as.character(seqnames(windows2))))
  peakView2 <- Views(peak.cov2[chr.idx2], as(windows2[strand(windows2)=="-"], "IntegerRangesList")[chr.idx2])## showMethods("coerce") can get all avaliable transforms
  
  peakView1 <- peakView1[lapply(peakView1, length)>0]
  peakView2 <- peakView2[lapply(peakView2, length)>0]
  diff_strandList <- list(peakView1,peakView2)
  diffStrandMatrixList <- lapply(diff_strandList, forEachStrand)
  return(diffStrandMatrixList)
}
aroundTss <- function(genomeFile, upstream = NULL,downstream = NULL, windowSize, methy) {
  k = 1
  ifelse(methy=="6mA", grlist <- list(m6AGrange,m6AGrange2), grlist <- list(m4CGrange,m4CGrange2))
  for(j in grlist) {
    geneMethylationCount <- countOverlaps(promoter,j)
    # geneMethylationCount <- countOverlaps(promoter, m6AGrange)
    promoter_trimed <- promoter[names(promoter) %in% names(geneMethylationCount[geneMethylationCount!=0])]
    promoterSeq <- getSeq(FaFile(genomeFile), promoter_trimed)  # get sequences by a GRange
    tt <- getAroundTssMatrix(promoter, j)
    tt[[2]] <- t(apply(tt[[2]],1,rev))
    tt1 <- do.call(rbind,tt)
    ss <- colSums(tt1)
    if(k==1){
      ifelse(methy=="6mA",countBase <- consensusMatrix(promoterSeq)[1,], countBase <- consensusMatrix(promoterSeq)[2,]) 
    }else if(k==2){
      ifelse(methy=="6mA",countBase <- consensusMatrix(promoterSeq)[4,], countBase <- consensusMatrix(promoterSeq)[3,])
    }
    ss <- ss / countBase
    ss[which(is.nan(ss))] <- 0   # because YDJ gene annotations have no UTR, so 503 positions are all ATG
    pos <- value <- NULL
    tagCount <- data.frame(pos = c(-499:1000), value = ss)
    winPos <- matrix(tagCount$pos, ncol = windowSize, byrow = TRUE)
    winValue <- matrix(tagCount$value, ncol = windowSize, byrow = TRUE)
    xx = rowMeans(winPos)
    yy = rowMeans(winValue) * 100
    if (k == 1) {
      plot(xx,yy, type="l",col="#6b318a", xlab = "Distance from TSS(bp)",ylim=c(0.40, 0.60),ylab = paste(methy,"6mA methylation occupancy (%)",sep = " "), lwd=2)
      legend("topright", c("Coding strand","Template strand"), lty = c(1,1,2), lwd = 2, 
             col = c("#6b318a", "#48bfd9"), bty = "n", cex = 0.8, xpd = TRUE)
    }else if(k == 2){
      lines(xx,yy,col="#48bfd9",lwd=2)
    }
    k=k+1
  }
}
aroundTss(genomeFile = "/Users/alexwang/0data/Morchella/ydj.fasta",windowSize = 50, methy = "6mA")

############### plot gene Body profile
ge <- genes(txdb)
start(ge) <- start(ge)-800; end(ge) <- end(ge) + 800  # no UTR in Morchella gene model prediction, so add +-300 as TSS/TES
ge <- subsetByOverlaps(ge, gs.gr, type = "within")
ge <- ge[start(ge) > 0]
getGeneBodyMatrix_profile <- function(window,GRangeObj){
  peak.cov1 <- coverage(GRangeObj) # coverage return Rle Object
  cov.len1 <- elementNROWS(peak.cov1)
  cov.width1 <- GRanges(seqnames = names(cov.len1), IRanges(start = rep(1,length(cov.len1)), end = cov.len1))
  windows1 <- subsetByOverlaps(window, cov.width1, type = "within")  # reconfirm all promoter regions are in genome
  chr.idx1 <- intersect(names(peak.cov1), unique(as.character(seqnames(windows1))))
  peakView1 <- Views(peak.cov1[chr.idx1], as(windows1[strand(windows1)=="+"], "IntegerRangesList")[chr.idx1])## showMethods("coerce") can get all avaliable transforms
  peakView2 <- Views(peak.cov1[chr.idx1], as(windows1[strand(windows1)=="-"], "IntegerRangesList")[chr.idx1])## showMethods("coerce") can get all avaliable transforms
  peakView1 <- peakView1[lapply(peakView1, length)>0]
  peakView2 <- peakView2[lapply(peakView2, length)>0]
  diff_strandList <- list(peakView1,peakView2)
  diffStrandMatrixList = list()
  for(y in 1:2){
    forEachView <- function(z,binNum=70){   # for each gene: split gene length to 50 bins and calculate sum methylation counts
      as.vector(z)
      aa <- floor(length(z)/50)
      ifelse(y==1, split.factor <- c(rep(-10:-1, each=50), rep(1:50,each=floor((length(z)-1000)/50)), rep(50,each=(length(z)-1000)%%50), rep(51:60, each=50)),
             split.factor <- c(rep(-10:-1, each=50), rep(1,each=(length(z)-1000)%%50), rep(1:50,each=floor((length(z)-1000)/50)), rep(51:60, each=50))
      )
      binCountSum <- sapply(split(z, split.factor), sum)
      return(binCountSum)
    }
    k=diff_strandList[[y]]
    library(BiocParallel)
    MulticoreParam = MulticoreParam(workers = 4)
    tagMatrixList <-bplapply(k, function(x) t(viewApply(x, forEachView)), BPPARAM = MulticoreParam)
    tagMatrix <- do.call("rbind", tagMatrixList)
    rownames(tagMatrix) <- 1:nrow(tagMatrix)
    tagMatrix <- tagMatrix[rowSums(tagMatrix) != 0,]   # remove no methylation genes
    diffStrandMatrixList[[y]] = tagMatrix
  }
  return(diffStrandMatrixList)
}

m6A_matrix1_igstr <- getGeneBodyMatrix_profile(ge, m6AGrange_igstr)
# m4C_matrix1_igstr <- getGeneBodyMatrix_profile(ge, m4CGrange_igstr)

geneMethylationCount <- countOverlaps(ge, m6AGrange)      # remove no methelation genes
ge_trimed <- ge[names(ge) %in% names(geneMethylationCount[geneMethylationCount!=0])]
ge_trimed <- ge_trimed[start(ge_trimed)>0]
geneBodySeq <- getSeq(FaFile("/Users/alexwang/0data/Morchella/ydj.fasta"), ge_trimed)
writeXStringSet(geneBodySeq, "/Users/alexwang/0data/Morchella/6mA/gebody.fasta")
system("awk 'BEGIN{i=1}/>/{$0=\">\"i;i++;print$0;next}{print$0}' /Users/alexwang/0data/Morchella/6mA/gebody.fasta > /Users/alexwang/0data/Morchella/6mA/gebody2.fasta")
system('python3 /Users/alexwang/0data/Morchella/6mA/countAC_inWindow.py m6A')

acrossGeneBody_profile <- function(genomeFile, upstream = NULL,downstream = NULL) {
  # for (i in 1:2) {
  i = 1
  # ifelse(i==1, 
  tt <- m6A_matrix1_igstr
  # , tt <- m4C_matrix1_igstr)
  # ifelse(i==1, 
  cc <- read.csv("/Users/alexwang/0data/Morchella/6mA/m6A_genebody",header = F)
  # , cc <- read.csv("/home/wanghm/wanghm_R/methylation/m4C_genebody",header = F))
  tt[[2]] <- t(apply(tt[[2]],1,rev))
  tt1 <- do.call(rbind,tt)
  ss <- colSums(tt1)  #  tt1 is a matrix with 70 cols, ss are 70 values, each value represent mean methy site counts in bins
  countBase <- cc$V1
  ss <- ss / countBase
  ss[which(is.nan(ss))] <- 0   # because YDJ gene annotations have no UTR, so 503 positions are all ATG
  pos <- c(seq(-475,0,50), seq(20,2000,40), seq(2025,2500,50))
  if(i==1){
    plot(pos,ss, type="l",col="#157a6e", xlab = "Genebody",ylab = "6mA methylation occupancy",lwd=2,xaxt = 'n')
    axis(1,at = c(0,2000), labels = c("TSS", "TTS"))
    legend("topright", c("6mA"), lty = c(1,1,2), lwd = 2, 
           col = c("#157a6e"), bty = "n", cex = 0.8, xpd = TRUE)
  }
  # else{
  #   lines(pos,ss,col="#e2711d",lwd=2)
  # }
}
# }
acrossGeneBody_profile(genomeFile = "/Users/alexwang/0data/Morchella/ydj.fasta")
##############################################
#### gene family ####
library(ggtree)
tree <- read.tree("/Users/alexwang/0data/Morchella/orthofinder/Species_Tree/SpeciesTree_rooted.txt")

tree <- 
  ggtree(tree)  + geom_tiplab(align = T,linetype = "dotted", size=3) + xlim(0,0.23)+ 
  geom_point(color='black', size=0.8) 

#############################################
#### Second metabolic gene family lineage search####
#!/bin/bash
## download ensembl species list 
# wget http://ftp.ensemblgenomes.org/vol1/pub/fungi/release-48/species_EnsemblFungi.txt
## extract taxonomy ID one per line
#awk '{gsub(";","\n",$8);print $8 > $1}' ../sm_protein_blastp.tbl
## extract species names and family name
#for i in `ls evm.model*`;do
## extract family name of ensembl fungi list
# taxonkit lineage -r <(cut -f 4 ~/0data/db/fungi_pep48/species_EnsemblFungi.txt) | awk -F "\t" '$3~/species/{ll=split($2,a,";");print a[ll-2]}'

#taxonkit lineage ${i} | awk -F ";" 'BEGIN{OFS="\t"}$0!~/unclassified/{split($NF,aa, " ");if(NF==13||NF==14){x=$11}else{if(NF==15||NF==16){x=$13}};print aa[1]" "aa[2],x}' >> species.list
#taxonkit lineage ${i} | awk -F ";" 'BEGIN{OFS="\t"}$0!~/unclassified/{split($NF,aa, " ");if(NF==13||NF==14){x=$11}else{if(NF==15||NF==16){x=$13}};print aa[1]" "aa[2],x}' | sort | uniq > ${i}.species
#done
## extract family name
#less species_list.txt | grep -v "incertae sedis" | awk -F "\t" '{print $2}' | sort | uniq | grep ".*eae" > sm_family.txt
## build topology tree
#perl ~/1.script/asm/Morchella/tax2tree.pl sm_family.txt

library(ggtree)
tree <- read.tree("/Users/alexwang/0data/Morchella/sm/scls/sm_taxid/sm_family.iTol.tree")
# tree <- midpoint.root(tree)

## color clade
tree <- groupClade(tree, .node = c("Agaricomycetes","Pezizales","leotiomyceta"))
p1 <- ggtree(tree, aes(color=group), branch.length = 'none', layout = "slanted",ladderize=F) + 
  geom_tiplab(fontface = "bold.italic") + scale_color_manual(values=c("black", "goldenrod3", "forestgreen","purple")) + 
  xlim(c(0, 10)) + theme(legend.position = "none") + geom_nodelab()
tree$node.label
## build 
matrix <- tree$tip.label
fh <-  list.files("/Users/alexwang/0data/Morchella/sm/scls/sm_taxid",pattern = 'evm.model.*species')
for (i in fh) {
  res = rep(0, length(tree$tip.label))
  d <- read.table(paste0("/Users/alexwang/0data/Morchella/sm/scls/sm_taxid/",i), fill=T, sep = "\t")
  tt <- table(d$V2)
  for (k in 1:length(tt)) {
    idx = which(tree$tip.label == names(tt)[k])
    res[idx] = tt[k]
  }
  matrix <- cbind(matrix, rev(res))
}
mat <- matrix[,2:ncol(matrix)]
mat <- apply(mat,2,as.numeric)
colnames(mat) = 1:13
SM_class <- data.frame(class = factor(c("terpene", "NRPS-like", "NRPS-like", "T1PKS", "NRPS-like", "terpene", "terpene", 
                                        "NRPS", "NRPS-like", "terpene", "NRPS-like", "NRPS-like", "terpene")))
rownames(SM_class) = colnames(mat)

p2 <- pheatmap::pheatmap(mat, cluster_rows = F, cluster_cols = F, display_numbers = T, 
                         number_format = "%s", border_color = "NA", 
                         color = colorRampPalette(c("white","red"))(15), legend = F, 
                         annotation_col = SM_class, show_colnames = F
                         # annotation_colors = list(class = RColorBrewer::brewer.pal(4, "Set3"))[1]
)
library(ggplotify)
p.heat <- as.ggplot(p2)
cowplot::plot_grid(p1, p.heat, ncol=2)
#### CAZymes class ####
AA=45
CBM=5
GH=130
PL=17
GT=58
CE=16
cazy.df <- readxl::read_xlsx("/Users/alexwang/0data/Morchella/0manuscript/Tables.xlsx", sheet = 2)
cazy.df$group = rep(c("AA", "CBM", "CE", "GH", "GT", "PL"), times=c(AA, CBM, CE, GH, GT, PL))

GH.count = sort(table(cazy.df[cazy.df$group == "GH",]$Family), decreasing = T)
GT.count = sort(table(cazy.df[cazy.df$group == "GT",]$Family), decreasing = T)
AA.count = sort(table(cazy.df[cazy.df$group == "AA",]$Family), decreasing = T)
PL.count = sort(table(cazy.df[cazy.df$group == "PL",]$Family), decreasing = T)
CE.count = sort(table(cazy.df[cazy.df$group == "CE",]$Family), decreasing = T)
CBM.count = sort(table(cazy.df[cazy.df$group == "CBM",]$Family), decreasing = T)

all.df <- rbind(data.frame(GH.count),
                data.frame(GT.count),
                data.frame(AA.count),
                data.frame(PL.count),
                data.frame(CE.count),
                data.frame(CBM.count))
write.table(all.df, file = "scls_cazy.csv", quote = F, row.names = F, sep = ",")

cc = c(colorRampPalette(c("#6247aa", "#dec9e9"))(length(GH.count)),
       colorRampPalette(c("#718355", "#e9f5db"))(length(GT.count)),
       colorRampPalette(c("#ff0a54", "#fae0e4"))(length(AA.count)),
       colorRampPalette(c("#ff7b00", "#ffea00"))(length(CBM.count)),
       colorRampPalette(c("#00b4d8", "#ade8f4"))(length(CE.count)),
       colorRampPalette(c("#b49286", "#d7bea8"))(length(PL.count)))
barplot(c(GH.count,
          GT.count,
          AA.count,
          CBM.count,
          CE.count,
          PL.count), 
        axes = T, width = 5,space=0, xlim=c(0,20), col=cc, cex.names = 0.2, ylab = "CAZyme family", xlab = "Number of CAZyme", horiz = T, border = F)





###############################################
#### species specific gene function ####

scls.go <- read.table("/Users/alexwang/0data/Morchella/go_annot.txt",header = T, sep = "\t", fill = T)
table(scls.go$TERM)
scls.sp <-read.table("/Users/alexwang/0data/Morchella/sextelata_comp/SCLS.sp.geneID")[,1]
# Morimp.sp <-read.table("/Users/alexwang/0data/Morchella/sextelata_comp/Morimp.sp.geneID")[,1]
scls.sp.go <- scls.go[scls.go$Gene %in% scls.sp, ]
table(scls.sp.go$TERM)
barplot(table(scls.sp.go$TERM))

library(VennDiagram)
morimp.orthogroup <- read.table("/Users/alexwang/0data/Morchella/sex_imp_orthofinder/Orthogroups/Morimp.orthogroups")[,1]
scls.orthogroup <- read.table("/Users/alexwang/0data/Morchella/sex_imp_orthofinder/Orthogroups/scls.orthogroups")[,1]
# venn.diagram(list(A=scls.orthogroup,B=morimp.orthogroup),filename="venn.tiff",lwd=1,lty=2,
#              col=c('red','green'),fill=c('red','green'),cat.col=c('red','green'),rotation.degree=90)

#### SCLS vs NZTD180501373 LS regions and SNPs ####
scls.genomeCov <- read.table("/Users/alexwang/0data/Morchella/colin/ydj.genomecov.txt")
scls.genomeCov$V1 <- gsub("Scaffold_", "", scls.genomeCov$V1)
scls.genomeCov.gr <- regioneR::toGRanges(scls.genomeCov[scls.genomeCov$V4 != 0, ])
## plot
gs <- read.table(paste0("/Users/alexwang/0data/Morchella/ydj.genome.txt"), header = T)
gs.gr <- GRanges(seqnames = gs$chr, ranges = IRanges(start=1, end = gs$size), col="red")
pp = getDefaultPlotParams(plot.type = 1)
pp$leftmargin = 0.05
pp$topmargin = 80
pp$bottommargin = 15
pp$ideogramheight = 20
pp$data1inmargin = 10
pp$data1outmargin = 0
kp = plotKaryotype(genome = gs.gr, plot.params = pp, cex=0.6)
kpPlotRegions(kp, data = scls.genomeCov.gr, r0 = 0, r1 = 0.6, col = "DarkCyan", avoid.overlapping = F)
## genes

scls.LS.gr <- regioneR::toGRanges(scls.genomeCov[scls.genomeCov$V4 == 0, ])
scls.LS.gr <- GRanges(seqnames = gsub("^","Scaffold_",seqnames(scls.LS.gr)), 
                      ranges = IRanges(start = start(scls.LS.gr), end = end(scls.LS.gr)))
LS.geneID <- subsetByOverlaps(gene, scls.LS.gr, type = "within")$gene_id
LS.gene.go <- scls.go[scls.go$Gene %in% LS.geneID, c("Gene", "GOID", "ONTOLOGY", "TERM")]
nrow(LS.gene.go)
sort(table(LS.gene.go$TERM))
er <- enricher(LS.geneID, TERM2GENE = term2gene)

expand.geneID <- read.table("/Users/alexwang/0data/Morchella/sex_imp_orthofinder/Orthogroups/scls.expand.geneID")[,1] # 268 gene families 
sp.geneID <- read.table("/Users/alexwang/0data/Morchella/sex_imp_orthofinder/Orthogroups/scls.sp.geneID")[,1] # 118 gene families
library(clusterProfiler)
er <- enricher(expand.geneID, TERM2GENE = term2gene, TERM2NAME = term2name)
dotplot(er)
# genes with unction 
# awk -F "\t" 'NF>11 && $4!~/MobiDBLite/{print}' scls_interpro.tsv > scls_interpro_withFunction.tsv
gene_withFunction.ID <- read.table("/Users/alexwang/0data/Morchella/with_funtions.geneID")[, 1]
length(intersect(gene_withFunction.ID, LS.geneID))

## SNPs
diff.df <- read.table("/Users/alexwang/0data/Morchella/colin/nztd_ydj.snps", skip = 4)
snp.df <- diff.df[diff.df$V2 != "." & diff.df$V3 != ".", ]
indel.df <- diff.df[diff.df$V2 == "." | diff.df$V3 == ".", ]

snp.gr <- GRanges(seqnames = snp.df$V12,
                  IRanges(start = snp.df$V4, end = snp.df$V4))
indel.gr <- GRanges(seqnames = indel.df$V12,
                  IRanges(start = indel.df$V4, end = indel.df$V4))
gene.gr <- genes(txdb)
length(unique(c(subsetByOverlaps(gene.gr, snp.gr)$gene_id, subsetByOverlaps(gene.gr, indel.gr)$gene_id)))






