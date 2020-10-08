library(circlize)
library(Rsamtools)
library(GenomicFeatures)

## prepare data
ydj.fai <- read.table("/Users/alexwang/0data/Morchella/ydj.fasta.fai")
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

chro <- ydj.fai$V1
starts = rep(0, length(canu.gs))
ends = canu.gs
genoCir <- data.frame(chr=chro, start=starts, end=ends)
genoCir$chr <- as.vector(genoCir[,1])
#### start plot######################################### only show fisrt 28 scaffolds of Morchella genome
circos.clear()
circos.par(start.degree = 87, track.height = 0.02, cell.padding = c(0,0,0,0), gap.degree=c(rep(1,27), 5))
circos.genomicInitialize(data = genoCir[1:28,],
                         sector.names = sub("Scaffold_","",chro[1:28]),
                         labels.cex = 0.5, track.height = 0.05, plotType = "labels")
### a: ideagram of 28 Scaffolds
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
    idx = which(chro[1:28] == sector.index)
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
cazy.geneID <- read.table("/Users/alexwang/0data/Morchella/cazy/ydj_cazy.geneID", header = T)[,1]
cazy.gr <- ge[ge$gene_id %in% cazy.geneID]
bed.cazy <- data.frame(chr=seqnames(cazy.gr), start=start(cazy.gr), end=end(cazy.gr))
bed.cazy$value <- 1
circos.genomicTrackPlotRegion(bed.cazy, track.height = 0.05,ylim = c(0,1), bg.border="gray",
                              panel.fun = function(region, value, ...){
                                circos.genomicRect(region, value, col = "darkorange", border = "darkorange", lwd = 1)
                              })

### e: SM clusters
# sm.geneID <- read.table("/Users/alexwang/0data/Morchella/sm/antismash/ydj.antismash.geneID")[,1]
# sm.gr <- ge[ge$gene_id %in% sm.geneID]
# bed.sm <- data.frame(chr=seqnames(sm.gr), start=start(sm.gr), end=end(sm.gr))
# bed.sm$value <- 1
# circos.genomicTrackPlotRegion(bed.sm, track.height = 0.05,ylim = c(0,1), bg.border="lightgray",
#                               panel.fun = function(region, value, ...){
#                                 circos.genomicRect(region, value, col = "darkgreen", border = "darkgreen", lwd = 1)
#                               })

# custom_scale <- function(v, mi, ma){
#   (v - min(v)) * (ma-mi) / (max(v)-min(v)) + mi
# }
### e: pacBio coverage barplot
pb.df <- read.table("/Users/alexwang/0data/Morchella/coverage/pb.regions.bed")
kk <- pb.df$V4
kk[which(kk %in% boxplot.stats(kk)$out)] <- mean(kk)
pb.df$V4 <- kk
colnames(pb.df) <- c("chr", "start", "end", "value")
circos.genomicTrackPlotRegion(pb.df, track.height=0.05, bg.border=NA, ylim = c(0, max(kk)),
                              panel.fun = function(region, value, ...) {
                                circos.genomicRect(region, value, ybottom = 0, ytop.column = 1,
                                                   col = "blue", border = NA, lwd = 0.5)
                              })
### f: plot GC content
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
### g: gene ortholog
# blastp.df <- read.table("/Users/alexwang/0data/Morchella/gene_duplication/self.blastp.filter.tab")
# ge.df <- data.frame(chr=seqnames(ge), start = start(ge), end = end(ge), id = ge$gene_id)
# ge.df <- ge.df[ge.df$chr %in% paste0("Scaffold_", 1:28),]
# 
# g1 <- blastp.df$V1
# g2 <- blastp.df$V2
# # reorder dataframe with specific order by g1
# r1 <- ge.df[match(g1, ge.df$id), c(1,2,3)]
# r2 <- ge.df[match(g2, ge.df$id), c(1,2,3)]
# bind.df <- cbind(r1,r2)
# bind.df <- na.omit(bind.df)
# rr1 <- bind.df[,c(1,2,3)]
# rr2 <- bind.df[,c(4,5,6)]
# circos.genomicLink(rr1, rr2, border = NA, col =  "pink", lwd = 3)

### g: LTR duplication
LTR_cluster <- read.table("/Users/alexwang/0data/Morchella/LTR/Cluster-pred-chrAll.285.0.match")[,c(2,6)]
LTR_cluster <- LTR_cluster[!duplicated(LTR_cluster),]
library(stringr)
r1 <- str_split_fixed(LTR_cluster$V2, "[-_:]", 5)
r2 <- str_split_fixed(LTR_cluster$V6, "[-_:]", 5)
rr1 <- data.frame(chr = paste0(r1[,1], "_", r1[,2]), start = as.numeric(r1[,3]), end = as.numeric(r1[,4]))
rr2 <- data.frame(chr = paste0(r2[,1], "_", r2[,2]), start = as.numeric(r2[,3]), end = as.numeric(r2[,4]))
### do not display LTR located on 29:42 Scaffolds
remove_row.idx <- unique(c(which(as.numeric(r1[,2]) > 28), which(as.numeric(r2[,2]) > 28)))
rr1 <- rr1[-remove_row.idx, ]
rr2 <- rr2[-remove_row.idx, ]
circos.genomicLink(rr1, rr2, border = NA, col =  "pink", lwd = 5)




