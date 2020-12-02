library(karyoploteR)
library(GenomicFeatures)
library(magrittr)
############ plot colin relationship between NZTD and SCLS
gs <- read.table(file = "/Users/alexwang/0data/pacbio/Morchella/colin/nztd.genome.txt", header = T)
gs <- gs[order(gs$size),]
gs.gr <- GRanges(seqnames = gs$chr, ranges = IRanges(start=1, end = gs$size), col="red")
## bundle regions
max.gap <- 5000
bundle1 <- read.table("/Users/alexwang/0data/pacbio/Morchella/colin/ydj.bundle_telo.txt", header = T)
bundle1 <- bundle1[abs(bundle1$s1-bundle1$e1)>max.gap & abs(bundle1$s2-bundle1$e2)>max.gap, ]
cc <- rainbow(20)[bundle1$ali]
x <- GRanges(seqnames = bundle1$ref, 
             ranges = IRanges(start = bundle1$s1, end = bundle1$e1),col = cc)
bundle.gr <- x
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

telo1 <- bundle1$telo1
telo2 <- bundle1$telo2
for (idx in 1:NROW(bundle.gr)){
  if (!is.na(telo1[idx])) {
    kp <- kpPoints(karyoplot = kp, chr = seqnames(bundle.gr[idx]), x = as.numeric(telo1[idx]), 
                   y=0.5, r0 = 0, r1 =0.8, ymax = 1, col="red", cex=0.4)
  }
  if (!is.na(telo2[idx])) {
    kp <- kpPoints(karyoplot = kp, chr = seqnames(bundle.gr[idx]), x = as.numeric(telo2[idx]), 
                   y=0.5, r0 =0 , r1 = 0.8, ymax = 1, col="red",cex = 0.4)
  }
}
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
txdb <- makeTxDbFromGFF("/Users/alexwang/0data/Morchella/ydj.gtf", format = "gtf")
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
#### start plot###############################################################
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







