library(ggtree)
tree <- read.tree("/Users/alexwang/0data/Morchella/sm/sm_taxid/sm_family.iTol.tree")
# tree <- midpoint.root(tree)

## color clade
tree <- groupClade(tree, .node = c("Agaricomycetes","Pezizales","leotiomyceta"))
p1 <- ggtree(tree, aes(color=group), branch.length = 'none', layout = "slanted",ladderize=F) + 
    geom_tiplab(fontface = "bold.italic") + scale_color_manual(values=c("black", "goldenrod3", "forestgreen","purple")) + 
    xlim(c(0, 10)) + theme(legend.position = "none") + geom_nodelab()
tree$node.label
## build 
matrix <- tree$tip.label
fh <-  list.files("/Users/alexwang/0data/Morchella/sm/sm_taxid",pattern = 'evm.model.*species')
for (i in fh) {
  res = rep(0, length(tree$tip.label))
  d <- read.table(paste0("/Users/alexwang/0data/Morchella/sm/sm_taxid/",i), fill=T, sep = "\t")
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


