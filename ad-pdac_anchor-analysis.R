setwd('C://Users/Lukas Vlahos/OneDrive/ac_lab/ad_pdac/')
library(Seurat)
library(PISCES)
library(ggplot2)
library(ComplexHeatmap)
library(reshape2)


### load general objects, marker sets; set colors
###############
sample.table <- read.table('sample-sheet.tsv', sep = '\t', header = TRUE)
sample.names <- c('KA003', 'KA005', 'KA006', 'KA007', 'KA008', 'KA009')
cluster.mappings <- read.csv('combined_analysis/cluster-mapping_20230124.csv')
## set colors
colnames(cluster.mappings) <- c('Cluster', 'Cell.Type', 'Sub.Type')
treatment.cols <- c('#000000', '#118002')
names(treatment.cols) <- c('Vehicle', '926')
cluster.cols <- group_colors(length(unique(cluster.mappings$Cluster)))
names(cluster.cols) <- sort(unique(cluster.mappings$Cluster))
ct.cols <- c('#107F80', '#66CCFE', '#FF0066', '#AA66FF', '#40007F')
names(ct.cols) <- c('Endothelial', 'Epithelial', 'Fibroblast', 'Lymphocyte', 'Myeloid')
sig.cols <- c('red', 'darkgrey', 'blue')
names(sig.cols) <- c('UP', 'NS', 'DOWN')
sig.shapes <- c(16, 15, 16)
names(sig.shapes) <- c('UP', 'NS', 'DOWN')
## marker sets
fibro.genes <- c('Gli1', 'Gli2', 'Gli3', 'Wif1', 'Ptch1')
wnt.receptors <- c('Lrp5', 'Lrp6', 'Ror1', 'Ror2', 'Ryk',
                   paste('Fzd', 1:10, sep = ''))
wnt.tfs <- c('Nfact1-4', 'Tcf7', 'Tcf7l1', 'Tcf7l2', 'Lef1', 'Ctnnb1')
endo.genes <- c('Kdr', 'Notch1')
wnt.genes <- c('Wnt1', 'Wnt2', 'Wnt2b', 'Wnt3', 'Wnt3A', 
               'Wnt4', 'Wnt5a', 'Wnt5b', 'Wnt6', 'Wnt7a', 
               'Wnt7b', 'Wnt8a', 'Wnt8b', 'Wnt9a', 'Wnt9b', 
               'Wnt10a', 'Wnt10b', 'Wnt11', 'Wnt16')
angio.ligands <- c('Vegfa', 'Vegfb', 'Vegfc', 'Pgf', 
                   'Thbs1', 'Thbs2', 'Thbs3', 'Thbs4', 
                   'Angpt1', 'Angpt2', 'Angptl2', 'Angptl3', 'Angpt4', 'Angptl6')
angio.receptors <- c('Flt1', 'Flt4',
                     'Kdr', 'CD47', 'CD36', 'Tek', 'Tie1')
###############

### tables of top proteins
###############
tf.regs <- readRDS('../ARACNe3_ljv/a3_regulator-list_tf-cotf.rds')
sig.regs <- readRDS('../ARACNe3_ljv/a3_regulator-list_sig.rds')
tf.list <- c(tf.regs$mus$tf$GN, tf.regs$mus$cotf$GN)
sig.list <- c(sig.regs$mus$sig$GN)
## make table for each
tf.table.list <- list()
sig.table.list <- list()
num.prots <- 100
for (ct in names(ct.cols)) {
  print(ct)
  narnea.obj <- readRDS(paste('anchor-combined_analysis/ct_group-narnea/', ct, '_narnea.rds', sep = ''))
  # get top sig prots
  pval.vec <- pnorm(abs(narnea.obj$NES[,1]), lower.tail = FALSE)
  pval.vec <- p.adjust(pval.vec, method = 'BH')
  sig.prots <- which(pval.vec < 0.05)
  sort.pes.vec <- sort(narnea.obj$PES[names(sig.prots), 1], decreasing = TRUE)
  # create tf data frame
  ct.tf.prots <- intersect(names(sort.pes.vec), tf.list)[1:num.prots]
  ct.tf.df <- data.frame('Protein' = ct.tf.prots,
                         'p.val' = pval.vec[ct.tf.prots],
                         'PES' = sort.pes.vec[ct.tf.prots])
  write.csv(ct.tf.df, file = paste('anchor-combined_analysis/anchor_analysis/ct_mr-lists/',
                                   ct, '_top-tfs_926-vs-vehicle.csv', sep = ''),
            quote = FALSE, row.names = FALSE)
  # create sig data frame
  ct.sig.prots <- intersect(names(sort.pes.vec), sig.list)[1:num.prots]
  ct.sig.df <- data.frame('Protein' = ct.sig.prots,
                          'p.val' = pval.vec[ct.sig.prots],
                          'PES' = sort.pes.vec[ct.sig.prots])
  write.csv(ct.sig.df, file = paste('anchor-combined_analysis/anchor_analysis/ct_mr-lists/',
                                    ct, '_top-sig_926-vs-vehicle.csv', sep = ''),
            quote = FALSE, row.names = FALSE)
}
## format to tables
tf.table <- Reduce(cbind, tf.table.list)
colnames(tf.table) <- names(tf.table.list)
sig.table <- Reduce(cbind, sig.table.list)
colnames(sig.table) <- names(sig.table.list)
## write as csvs
write.csv(tf.table, file = 'anchor-combined_analysis/anchor_analysis/ct_top-tf-prots_926-vs-vehicle.csv',
          quote = FALSE, row.names = FALSE, col.names = TRUE)
write.csv(sig.table, file = 'anchor-combined_analysis/anchor_analysis/ct_top-signaling-prots_926-vs-vehicle.csv',
          quote = FALSE, row.names = FALSE)
###############

### tables of bot proteins
###############
tf.regs <- readRDS('../ARACNe3_ljv/a3_regulator-list_tf-cotf.rds')
sig.regs <- readRDS('../ARACNe3_ljv/a3_regulator-list_sig.rds')
tf.list <- c(tf.regs$mus$tf$GN, tf.regs$mus$cotf$GN)
sig.list <- c(sig.regs$mus$sig$GN)
## make table for each
tf.table.list <- list()
sig.table.list <- list()
num.prots <- 100
for (ct in names(ct.cols)) {
  print(ct)
  narnea.obj <- readRDS(paste('anchor-combined_analysis/ct_group-narnea/', ct, '_narnea.rds', sep = ''))
  # get top sig prots
  pval.vec <- pnorm(abs(narnea.obj$NES[,1]), lower.tail = FALSE)
  pval.vec <- p.adjust(pval.vec, method = 'BH')
  sig.prots <- which(pval.vec < 0.05)
  sort.pes.vec <- sort(narnea.obj$PES[names(sig.prots), 1], decreasing = FALSE )
  # create tf data frame
  ct.tf.prots <- intersect(names(sort.pes.vec), tf.list)[1:num.prots]
  ct.tf.df <- data.frame('Protein' = ct.tf.prots,
                         'p.val' = pval.vec[ct.tf.prots],
                         'PES' = sort.pes.vec[ct.tf.prots])
  write.csv(ct.tf.df, file = paste('anchor-combined_analysis/anchor_analysis/ct_mr-lists/',
                                   ct, '_bot-tfs_926-vs-vehicle.csv', sep = ''),
            quote = FALSE, row.names = FALSE)
  # create sig data frame
  ct.sig.prots <- intersect(names(sort.pes.vec), sig.list)[1:num.prots]
  ct.sig.df <- data.frame('Protein' = ct.sig.prots,
                          'p.val' = pval.vec[ct.sig.prots],
                          'PES' = sort.pes.vec[ct.sig.prots])
  write.csv(ct.sig.df, file = paste('anchor-combined_analysis/anchor_analysis/ct_mr-lists/',
                                    ct, '_bot-sig_926-vs-vehicle.csv', sep = ''),
            quote = FALSE, row.names = FALSE)
}
###############


### 1 - scatter plot colored by treatment
###############
seurat.obj <- readRDS('anchor-combined_analysis/anchored_seurat-obj.rds')
cell.names <- colnames(seurat.obj)
cell.samples <- sapply(cell.names, function(x) {strsplit(x, '\\.')[[1]][1]} )
cell.treatment <- sample.table$Group[match(cell.samples, sample.table$Sample)]
## make plot df
plot.df <- data.frame('UMAP1' = seurat.obj@reductions$umap@cell.embeddings[,1],
                      'UMAP2' = seurat.obj@reductions$umap@cell.embeddings[,2],
                      'Treatment' = cell.treatment)
cluster.umap.centers <- lapply(sort(unique(seurat.obj$seurat_clusters)), function(x) {
  clust.samps <- which(seurat.obj$seurat_clusters == x)
  mean.vals <- colMeans(seurat.obj@reductions$umap@cell.embeddings[clust.samps, 1:2])
  return(mean.vals)
})
cluster.umap.centers <- Reduce(rbind, cluster.umap.centers)
cluster.umap.centers <- as.data.frame(cluster.umap.centers)
cluster.umap.centers$Cluster <- sort(unique(seurat.obj$seurat_clusters))
## make plot
ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = Treatment), size = 0.5) +
  ggtitle('Combined Analysis - Treatment') +
  scale_color_manual(values = treatment.cols) + 
  guides(colour = guide_legend(override.aes = list(size = 3)))
ggsave('anchor-combined_analysis/anchor_analysis/fig-1_gexp-umap-treatment.pdf',
       height = 8, width = 10, units = 'in', dpi = 300)
###############

### 2a/b - scatter plot w/ cluster and cell type
###############
seurat.obj <- readRDS('anchor-combined_analysis/anchored_seurat-obj.rds')
cell.names <- colnames(seurat.obj)
cell.samples <- sapply(cell.names, function(x) {strsplit(x, '\\.')[[1]][1]} )
cell.treatment <- sample.table$Group[match(cell.samples, sample.table$Sample)]
clust.vec <- seurat.obj$seurat_clusters
ct.vec <- plyr::mapvalues(clust.vec, from = cluster.mappings$Cluster, to = cluster.mappings$Cell.Type)
ct.vec <- plyr::mapvalues(ct.vec, from = sort(as.character(unique(ct.vec))),
                          to = sort(names(ct.cols)))
## make plot df
plot.df <- data.frame('UMAP1' = seurat.obj@reductions$umap@cell.embeddings[,1],
                      'UMAP2' = seurat.obj@reductions$umap@cell.embeddings[,2],
                      'Cluster' = as.factor(clust.vec),
                      'Cell.Type' = as.factor(ct.vec))
# cluster centers
cluster.umap.centers <- lapply(sort(unique(seurat.obj$seurat_clusters)), function(x) {
  clust.samps <- which(seurat.obj$seurat_clusters == x)
  mean.vals <- colMeans(seurat.obj@reductions$umap@cell.embeddings[clust.samps, 1:2])
  return(mean.vals)
})
cluster.umap.centers <- Reduce(rbind, cluster.umap.centers)
cluster.umap.centers <- as.data.frame(cluster.umap.centers)
cluster.umap.centers$Cluster <- sort(unique(seurat.obj$seurat_clusters))
# ct centers
ct.umap.centers <- lapply(sort(unique(ct.vec)), function(x) {
  ct.samps <- which(ct.vec == x)
  mean.vals <- colMeans(seurat.obj@reductions$umap@cell.embeddings[ct.samps, 1:2])
  return(mean.vals)
})
ct.umap.centers <- Reduce(rbind, ct.umap.centers)
ct.umap.centers <- as.data.frame(ct.umap.centers)
ct.umap.centers$ct <- sort(unique(ct.vec))
## make plot - cluster
ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = Cluster), size = 0.5) +
  ggtitle('Combined Analysis - Cluster') +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  scale_color_manual(values = cluster.cols) + 
  annotate("text", x = cluster.umap.centers$UMAP_1, y = cluster.umap.centers$UMAP_2,
           label = cluster.umap.centers$Cluster)
ggsave('anchor-combined_analysis/anchor_analysis/fig-2-a_gexp-umap-cluster.pdf',
       height = 8, width = 10, units = 'in', dpi = 300)
## make plot - cell type
ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = Cell.Type), size = 0.5) +
  ggtitle('Combined Analysis - Cell Type') +
  guides(colour = guide_legend(override.aes = list(size = 3), title = 'Cell Type')) +
  scale_color_manual(values = ct.cols) + 
  annotate("text", x = ct.umap.centers$UMAP_1, y = ct.umap.centers$UMAP_2,
           label = ct.umap.centers$ct, size = 8)
ggsave('anchor-combined_analysis/anchor_analysis/fig-2-b_gexp-umap-ct.pdf',
       height = 8, width = 10, units = 'in', dpi = 300)
###############

### 3a/b - bar graph of ct frequencies for treatment / sample
###############
seurat.obj <- readRDS('anchor-combined_analysis/anchored_seurat-obj.rds')
cell.names <- colnames(seurat.obj)
cell.samples <- sapply(cell.names, function(x) {strsplit(x, '\\.')[[1]][1]} )
cell.treatment <- sample.table$Group[match(cell.samples, sample.table$Sample)]
names(cell.treatment) <- names(cell.samples)
clust.vec <- seurat.obj$seurat_clusters
ct.vec <- plyr::mapvalues(clust.vec, from = cluster.mappings$Cluster, to = cluster.mappings$Cell.Type)
ct.vec <- plyr::mapvalues(ct.vec, from = sort(as.character(unique(ct.vec))),
                          to = sort(names(ct.cols)))
## treatment bargraph
conf.table <- table(cell.treatment, ct.vec)
conf.table <- conf.table / rowSums(conf.table)
plot.df <- melt(conf.table)
colnames(plot.df) <- c('Treatment', 'Cell.Type', 'Percentage')
ggplot(plot.df, aes(x = Treatment, group = Cell.Type)) + 
  geom_bar(aes(y = Percentage, fill = Cell.Type), stat = "identity", position = position_dodge()) +
  labs(y = "Cell Type Percentage", x = "Treatment", title = "Cell Type Frequency by Treatment") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = ct.cols) +
  geom_text(aes(label = scales::percent(Percentage), y = Percentage, group = Cell.Type), 
            position = position_dodge(width = 0.9), vjust = -.5)
ggsave('anchor-combined_analysis/anchor_analysis/fig-3-a_ct-treatment-bargraph.pdf',
       height = 6, width = 10, units = 'in', dpi = 300)
# stackged bargraph
sort.ct.names <- sort(as.character(unique(plot.df$Cell.Type)))
plot.df$Cell.Type <- factor(plot.df$Cell.Type, levels = sort.ct.names)
rev.sort.ct.names <- sort(as.character(unique(plot.df$Cell.Type)), decreasing = TRUE)
pos.vec <- apply(plot.df, 1, function(x) {
  group.mat <- plot.df[which(plot.df$Treatment == x[1]),]
  group.mat <- group.mat[order(as.character(group.mat$Cell.Type), decreasing = TRUE),]
  group.mat <- group.mat[1:which(rev.sort.ct.names == x[2]),,drop = FALSE]
  x.pos <- sum(as.numeric(group.mat$Percentage)) - (0.5) * as.numeric(x[3])
  print(x.pos)
  return(x.pos)
})
plot.df$text.pos <- pos.vec
ggplot(plot.df, aes(x = Treatment, group = Cell.Type)) + 
  geom_bar(aes(y = Percentage, fill = Cell.Type), stat = "identity") +
  labs(y = "Cell Type Percentage", x = "Treatment", title = "Cell Type Frequency by Treatment") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = ct.cols) +
  geom_text(aes(label = scales::percent(Percentage), y = text.pos, group = Cell.Type))
ggsave('anchor-combined_analysis/anchor_analysis/fig-3-a_ct-treatment-bargraph_stacked.pdf',
       height = 6, width = 10, units = 'in', dpi = 300)
## treatment bargraph
conf.table <- table(cell.samples, ct.vec)
conf.table <- conf.table / rowSums(conf.table)
plot.df <- melt(conf.table)
colnames(plot.df) <- c('Sample', 'Cell.Type', 'Percentage')
ggplot(plot.df, aes(x = Sample, group = Cell.Type)) + 
  geom_bar(aes(y = Percentage, fill = Cell.Type), stat = "identity", position = position_dodge()) +
  labs(y = "Cell Type Percentage", x = "Sample", title = "Cell Type Frequency by Sample") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = ct.cols) +
  geom_text(aes(label = scales::percent(round(Percentage, digits = 3)), y = Percentage, group = Cell.Type), 
            position = position_dodge(width = 0.9), vjust = -.5)
ggsave('anchor-combined_analysis/anchor_analysis/fig-3-b_ct-sample-bargraph.pdf',
       height = 6, width = 10, units = 'in', dpi = 300)
# stackged bargraph
sort.ct.names <- sort(as.character(unique(plot.df$Cell.Type)))
plot.df$Cell.Type <- factor(plot.df$Cell.Type, levels = sort.ct.names)
rev.sort.ct.names <- sort(as.character(unique(plot.df$Cell.Type)), decreasing = TRUE)
pos.vec <- apply(plot.df, 1, function(x) {
  group.mat <- plot.df[which(plot.df$Sample == x[1]),]
  group.mat <- group.mat[order(as.character(group.mat$Cell.Type), decreasing = TRUE),]
  group.mat <- group.mat[1:which(rev.sort.ct.names == x[2]),,drop = FALSE]
  x.pos <- sum(as.numeric(group.mat$Percentage)) - (0.5) * as.numeric(x[3])
  print(x.pos)
  return(x.pos)
})
plot.df$text.pos <- pos.vec
ggplot(plot.df, aes(x = Sample, group = Cell.Type)) + 
  geom_bar(aes(y = Percentage, fill = Cell.Type), stat = "identity") +
  labs(y = "Cell Type Percentage", x = "Sample", title = "Cell Type Frequency by Sample") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = ct.cols) +
  geom_text(aes(label = scales::percent(round(Percentage, digits = 3)), y = text.pos, group = Cell.Type))
ggsave('anchor-combined_analysis/anchor_analysis/fig-3-b_ct-sample-bargraph_stacked.pdf',
       height = 6, width = 10, units = 'in', dpi = 300)
###############

### 4 - dot plot for differential expression analysis of fibroblast markers
###############
seurat.obj <- readRDS('anchor-combined_analysis/anchored_seurat-obj.rds')
cell.names <- colnames(seurat.obj)
cell.samples <- sapply(cell.names, function(x) {strsplit(x, '\\.')[[1]][1]} )
cell.treatment <- sample.table$Group[match(cell.samples, sample.table$Sample)]
names(cell.treatment) <- names(cell.samples)
clust.vec <- seurat.obj$seurat_clusters
ct.vec <- plyr::mapvalues(clust.vec, from = cluster.mappings$Cluster, to = cluster.mappings$Cell.Type)
rm(seurat.obj)
## load count matrix
cpm.mat <- readRDS('anchor-combined_analysis/combined-gexp_pad-filt-counts_cpm.rds')
## dif gene exp analysis for each cell type
dif.exp.mats <- list()
for (ct in unique(ct.vec)) {
  print(ct)
  # get group labels
  ct.cells <- names(ct.vec)[which(ct.vec == ct)]
  ct.926.cells <- intersect(ct.cells, names(cell.treatment)[which(cell.treatment == '926')])
  ct.vehicle.cells <- intersect(ct.cells, names(cell.treatment)[which(cell.treatment == 'Vehicle')])
  # perform analysis
  dif.exp <- sapply(fibro.genes, function(x) {
    test.vec <- cpm.mat[x, ct.926.cells]
    ref.vec <- cpm.mat[x, ct.vehicle.cells]
    test.obj <- wilcox.test(test.vec, ref.vec, alternative = "two.sided")
    p.val <- test.obj$p.val
    rbs.cor <- 2 * test.obj$statistic / (length(ct.926.cells) * length(ct.vehicle.cells)) - 1
    return(c(x, p.val, rbs.cor))
  })
  # reformat to matrix
  rownames(dif.exp) <- c('gene', 'pval', 'rbsc')
  dif.exp.mat <- as.data.frame(t(dif.exp))
  dif.exp.mat$ct <- ct
  dif.exp.mat$pval <- as.numeric(dif.exp.mat$pval)
  dif.exp.mat$rbsc <- as.numeric(dif.exp.mat$rbsc)
  dif.exp.mat$logp <- log(dif.exp.mat$pval)
  sig.vec <- rep('NS', nrow(dif.exp.mat))
  sig.vec[which(dif.exp.mat$pval < 0.05)] <- 'DOWN'
  dif.exp.mat$sig <- sig.vec
  # calculate NES
  nes.vec <- qnorm(dif.exp.mat$logp, log.p = TRUE, lower.tail = FALSE) * sign(dif.exp.mat$rbsc)
  dif.exp.mat$nes <- nes.vec
  # add to list
  dif.exp.mats[[ct]] <- dif.exp.mat
}
## pooled differential expression
pooled.926.cells <- names(cell.treatment)[which(cell.treatment == '926')]
pooled.vehicle.cells <- names(cell.treatment)[which(cell.treatment == 'Vehicle')]
# perform analysis
dif.exp <- sapply(fibro.genes, function(x) {
  test.vec <- cpm.mat[x, pooled.926.cells]
  ref.vec <- cpm.mat[x, pooled.vehicle.cells]
  test.obj <- wilcox.test(test.vec, ref.vec, alternative = "two.sided")
  p.val <- test.obj$p.val
  rbs.cor <- 2 * test.obj$statistic / (length(pooled.926.cells) * length(pooled.vehicle.cells)) - 1
  return(c(x, p.val, rbs.cor))
})
# reformat to matrix
rownames(dif.exp) <- c('gene', 'pval', 'rbsc')
dif.exp.mat <- as.data.frame(t(dif.exp))
dif.exp.mat$ct <- 'Pseudobulk'
dif.exp.mat$pval <- as.numeric(dif.exp.mat$pval)
dif.exp.mat$rbsc <- as.numeric(dif.exp.mat$rbsc)
dif.exp.mat$logp <- log(dif.exp.mat$pval)
sig.vec <- rep('NS', nrow(dif.exp.mat))
sig.vec[which(dif.exp.mat$pval < 0.05)] <- 'DOWN'
dif.exp.mat$sig <- sig.vec
# calculate NES
nes.vec <- qnorm(dif.exp.mat$logp, log.p = TRUE, lower.tail = FALSE) * sign(dif.exp.mat$rbsc)
dif.exp.mat$nes <- nes.vec
dif.exp.mats[['pseudobulk']] <- dif.exp.mat
saveRDS(dif.exp.mats, file = 'anchor-combined_analysis/fibro-genes_dif-gexp_cpm-mat.rds')
## make fibroblast plot
fibro.mat <- dif.exp.mats$Fibroblast
ggplot(fibro.mat, aes(gene, ct)) + geom_point(aes(color = nes, size = abs(rbsc))) +
  labs(x = 'Gene', y = 'Cell Type', title = 'Fibroblast Gene Differential Expression') +
  #scale_color_manual(values = sig.cols) + 
  scale_color_gradient2(low = 'blue', mid = 'grey', high = 'red', midpoint = 0) +
  scale_size_continuous(range = c(3, 10), limits = c(0, 0.2)) + 
  guides(size = guide_legend(title="|PES|", order = 2),
         color = guide_colourbar(title = "NES", order = 1))
ggsave('anchor-combined_analysis/anchor_analysis/fig-4_fibroblast-dif-exp-dot-plot.pdf',
       height = 4, width = 10, units = 'in', dpi = 300)
## all cell types (shaped)
plot.mat <- Reduce(rbind, dif.exp.mats)
ggplot(plot.mat, aes(gene, ct)) + geom_point(aes(color = nes, size = abs(rbsc), shape = sig)) +
  labs(x = 'Gene', y = 'Cell Type', title = 'Fibroblast Gene Differential Expression - All Cell Types') +
  #scale_color_manual(values = sig.cols) + 
  scale_color_gradient2(low = 'blue', mid = 'grey', high = 'red', midpoint = 0) +
  scale_size_continuous(range = c(3, 10), limits = c(0, 0.2)) + 
  scale_shape_manual(values = sig.shapes) + 
  guides(size = guide_legend(title="|PES|", order = 2),
         color = guide_colourbar(title = "NES", order = 1),
         shape = guide_legend(title = 'Significance', override.aes = list(size = 3))) +
  theme(legend.box = 'horizontal')
ggsave('anchor-combined_analysis/anchor_analysis/fig-4s_fibro-genes-all-cells-dif-exp-dot-plot_shaped.pdf',
       height = 4, width = 10, units = 'in', dpi = 300)
## all cell types (border)
plot.mat <- Reduce(rbind, dif.exp.mats)
ggplot(plot.mat, aes(gene, ct)) + geom_point(aes(fill = nes, size = abs(rbsc), colour = sig), pch = 21) +
  labs(x = 'Gene', y = 'Cell Type', title = 'Fibroblast Gene Differential Expression - All Cell Types') +
  #scale_color_manual(values = sig.cols) + 
  scale_fill_gradient2(low = 'blue', mid = 'grey', high = 'red', midpoint = 0) +
  scale_size_continuous(range = c(3, 10), limits = c(0, 0.2)) + 
  scale_color_manual(values = sig.cols) + 
  guides(size = guide_legend(title="|PES|", order = 2, override.aes = list(fill = 'darkgrey')),
         fill = guide_colourbar(title = "NES", order = 1),
         colour = guide_legend(title = 'Significance', override.aes = list(size = 3))) +
  theme(legend.box = 'horizontal')
ggsave('anchor-combined_analysis/anchor_analysis/fig-4s_fibro-genes-all-cells-dif-exp-dot-plot_border.pdf',
       height = 4, width = 10, units = 'in', dpi = 300)
## all cell types (binned)
plot.mat <- Reduce(rbind, dif.exp.mats)
sig.plot.mat <- plot.mat[which(plot.mat$sig != 'NS'),]
ns.plot.mat <- plot.mat[which(plot.mat$sig == 'NS'),]
ggplot(sig.plot.mat, aes(gene, ct)) + geom_point(aes(color = nes, size = abs(rbsc))) +
  labs(x = 'Gene', y = 'Cell Type', title = 'Fibroblast Gene Differential Expression - All Cell Types') +
  #scale_color_manual(values = sig.cols) + 
  scale_color_gradient2(low = 'blue', mid = 'grey', high = 'red', midpoint = 0) +
  scale_size_continuous(range = c(3, 10), limits = c(0, 0.2)) + 
  guides(size = guide_legend(title="|PES|", order = 2),
         color = guide_colourbar(title = "NES", order = 1)) +
  theme(legend.box = 'horizontal') + theme_bw() +
  geom_point(data = ns.plot.mat, aes(gene, ct, color = nes, size = abs(rbsc)), color = 'black')
ggsave('anchor-combined_analysis/anchor_analysis/fig-4s_fibro-genes-all-cells-dif-exp-dot-plot_binned.pdf',
       height = 4, width = 10, units = 'in', dpi = 300)
###############

### 5 - dot plot for differential activity analysis of fibroblast markers
###############
narnea.obj <- readRDS('anchor-combined_analysis/ct_group-narnea/Fibroblast_narnea.rds')
## prepare plot data
use.prots <- intersect(rownames(narnea.obj$NES), fibro.genes)
nes.vec <- narnea.obj$NES[use.prots,]
pes.vec <- narnea.obj$PES[use.prots,]
# convert to p-values
pval.vec <- pnorm(abs(nes.vec), lower.tail = FALSE)
sig.vec <- rep('NS', length(use.prots)); names(sig.vec) <- use.prots
sig.vec[intersect(which(pval.vec < 0.05), which(pes.vec > 0))] <- 'UP'
sig.vec[intersect(which(pval.vec < 0.05), which(pes.vec < 0))] <- 'DOWN'
# make plot data
plot.df <- data.frame('sig' = sig.vec,
                      'pes' = abs(pes.vec),
                      'nes' = nes.vec,
                      'gene' = use.prots,
                      'ct' = 'Fibroblast')
# make plot
sig.plot.mat <- plot.df[which(plot.df$sig != 'NS'),]
ns.plot.mat <- plot.df[which(plot.df$sig == 'NS'),]
ggplot(sig.plot.mat, aes(gene, ct)) + geom_point(aes(color = nes, size = abs(pes))) +
  labs(x = 'Protein', y = 'Cell Type', title = 'Fibroblast Protein Differential Activity') +
  #scale_color_manual(values = sig.cols) + 
  scale_color_gradient2(low = 'blue', mid = 'grey', high = 'red', midpoint = 0) +
  scale_size_continuous(range = c(3, 10), limits = c(0, 0.2)) + 
  guides(size = guide_legend(title="|PES|", order = 2),
         color = guide_colourbar(title = "NES", order = 1)) +
  theme_bw() +
  geom_point(data = ns.plot.mat, aes(gene, ct, color = nes, size = abs(pes)), color = 'black')
ggsave('anchor-combined_analysis/anchor_analysis/fig-5_fibroblast-dif-activity-dot-plot.pdf',
       height = 4, width = 10, units = 'in', dpi = 300)
###############

### 5s - dot plots for differntial activity analysis of fibro markers in all cell types
###############
ct.plot.dfs <- list()
for (ct in names(ct.cols)) {
  print(ct)
  narnea.obj <- readRDS(paste('anchor-combined_analysis/ct_group-narnea/', ct, '_narnea.rds', sep = ''))
  # get vectors
  use.prots <- intersect(rownames(narnea.obj$NES), fibro.genes)
  nes.vec <- narnea.obj$NES[use.prots,]
  pes.vec <- narnea.obj$PES[use.prots,]
  # convert to significance
  pval.vec <- pnorm(abs(nes.vec), lower.tail = FALSE)
  sig.vec <- rep('NS', length(use.prots)); names(sig.vec) <- use.prots
  sig.vec[intersect(which(pval.vec < 0.05), which(pes.vec > 0))] <- 'UP'
  sig.vec[intersect(which(pval.vec < 0.05), which(pes.vec < 0))] <- 'DOWN'
  # make plot data
  plot.df <- data.frame('sig' = sig.vec,
                        'pes' = abs(pes.vec),
                        'nes' = nes.vec,
                        'gene' = use.prots,
                        'ct' = ct)
  ct.plot.dfs[[ct]] <- plot.df
}
## make plot
plot.df <- Reduce(rbind, ct.plot.dfs)
sig.plot.mat <- plot.df[which(plot.df$sig != 'NS'),]
ns.plot.mat <- plot.df[which(plot.df$sig == 'NS'),]
ggplot(plot.df, aes(gene, ct)) + geom_point(aes(color = nes, size = abs(pes))) +
  labs(x = 'Protein', y = 'Cell Type', title = 'Fibroblast Protein Differential Activity - All Cell Types') +
  #scale_color_manual(values = sig.cols) + 
  scale_color_gradient2(low = 'blue', mid = 'grey', high = 'red', midpoint = 0) +
  scale_size_continuous(range = c(3, 10), limits = c(0, 0.2)) + 
  guides(size = guide_legend(title="|PES|", order = 2),
         color = guide_colourbar(title = "NES", order = 1)) +
  theme_bw() +
  geom_point(data = ns.plot.mat, aes(gene, ct, color = nes, size = abs(pes)), color = 'black')
ggsave('anchor-combined_analysis/anchor_analysis/fig-5s_fibroblast-dif-activity-dot-plot.pdf',
       height = 4,  width = 10, units = 'in', dpi = 300)
##############

### 6 + 6s - dot plot for differential activity analysis of wnt receptors
###############
ct.receptor.dfs <- list()
for (ct in sort(unique(ct.vec))) {
  cat(paste(ct, '\n'))
  narnea.res <- readRDS(paste('anchor-combined_analysis/ct_group-narnea/', ct, '_narnea.rds', sep = ''))
  receptor.intersect <- intersect(rownames(narnea.res$NES), wnt.receptors)
  receptor.nes <- narnea.res$NES[receptor.intersect,]
  receptor.pes <- narnea.res$PES[receptor.intersect,]
  ct.receptor.df <- data.frame('Cell.Type' = ct,
                               'Protein' = receptor.intersect,
                               'PES' = receptor.pes,
                               'NES' = receptor.nes)
  ct.receptor.dfs[[ct]] <- ct.receptor.df
}
dot.plot.df <- Reduce(rbind, ct.receptor.dfs)
## reformat dataframe
dot.plot.df <- as.data.frame(dot.plot.df)
dot.plot.df$Cell.Type <- as.factor(dot.plot.df$Cell.Type)
dot.plot.df$Protein <- as.factor(dot.plot.df$Protein)
dot.plot.df$PES <- abs(dot.plot.df$PES)
## calculate significance
pval.vec <- pnorm(abs(dot.plot.df$NES), lower.tail = FALSE)
pval.vec <- p.adjust(pval.vec, method = 'BH')
names(pval.vec) <- rownames(dot.plot.df)
sig.vec <- rep('NS', length(pval.vec))
names(sig.vec) <- rownames(dot.plot.df)
sig.vec[intersect(names(which(pval.vec < 0.05)), rownames(dot.plot.df)[which(dot.plot.df$NES > 0)])] <- 'UP'
sig.vec[intersect(names(which(pval.vec < 0.05)), rownames(dot.plot.df)[which(dot.plot.df$NES < 0)])] <- 'DOWN'
dot.plot.df$SIG <- as.factor(sig.vec)
## make plot - supplement
sig.plot.mat <- dot.plot.df[which(dot.plot.df$SIG != 'NS'),]
ns.plot.mat <- dot.plot.df[which(dot.plot.df$SIG == 'NS'),]
ggplot(sig.plot.mat, aes(x = Protein, y = Cell.Type)) +
  geom_point(aes(size = PES, color = NES)) +
  labs(title = 'WNT Receptor Activity in 926 vs Vehicle', y = 'Cell Type')  +
  #scale_color_manual(values = sig.cols) + 
  scale_color_gradient2(low = 'blue', mid = 'grey', high = 'red', midpoint = 0) +
  scale_size_continuous(range = c(3, 10), limits = c(0, 0.2)) + 
  guides(size = guide_legend(title="|PES|", order = 2),
         color = guide_colourbar(title = "NES", order = 1)) +
  theme_bw() +
  geom_point(data = ns.plot.mat, aes(Protein, Cell.Type, color = NES, size = abs(PES)), color = 'black')
ggsave('anchor-combined_analysis/anchor_analysis/fig-6s_wnt-receptor-dif-activity-all.pdf',
       height = 4, width = 10, units = 'in', dpi = 300)
## susbet
main.subset <- c('Ryk', 'Ror1', 'Ror2', 'Lrp5', 'Lrp6', 'Fzd1', 'Fzd2', 'Fzd3', 'Fzd4')
main.plot.df <- dot.plot.df[which(dot.plot.df$Protein %in% main.subset),]
sig.plot.mat <- main.plot.df[which(main.plot.df$SIG != 'NS'),]
ns.plot.mat <- main.plot.df[which(main.plot.df$SIG == 'NS'),]
ggplot(sig.plot.mat, aes(x = Protein, y = Cell.Type)) +
  geom_point(aes(size = PES, color = NES)) + 
  labs(title = 'WNT Receptor Activity in 926 vs Vehicle', y = 'Cell Type')  +
  #scale_color_manual(values = sig.cols) + 
  scale_color_gradient2(low = 'blue', mid = 'grey', high = 'red', midpoint = 0) +
  scale_size_continuous(range = c(3, 10), limits = c(0, 0.2)) + 
  guides(size = guide_legend(title="|PES|", order = 2),
         color = guide_colourbar(title = "NES", order = 1)) +
  theme_bw() +
  geom_point(data = ns.plot.mat, aes(Protein, Cell.Type, color = NES, size = abs(PES)), color = 'black')
ggsave('anchor-combined_analysis/anchor_analysis/fig-6_wnt-receptor-dif-activity.pdf',
       height = 4, width = 10, units = 'in', dpi = 300)
###############

### 7 - dot plot for differential activity analysis of wnt tfs
###############
ct.tf.dfs <- list()
for (ct in sort(unique(ct.vec))) {
  cat(paste(ct, '\n'))
  narnea.res <- readRDS(paste('anchor-combined_analysis/ct_group-narnea/', ct, '_narnea.rds', sep = ''))
  tf.intersect <- intersect(rownames(narnea.res$NES), wnt.tfs)
  tf.nes <- narnea.res$NES[tf.intersect,]
  tf.pes <- narnea.res$PES[tf.intersect,]
  ct.tf.df <- data.frame('Cell.Type' = ct,
                         'Protein' = tf.intersect,
                         'PES' = tf.pes,
                         'NES' = tf.nes)
  ct.tf.dfs[[ct]] <- ct.tf.df
}
dot.plot.df <- Reduce(rbind, ct.tf.dfs)
## reformat dataframe
dot.plot.df <- as.data.frame(dot.plot.df)
dot.plot.df$Cell.Type <- as.factor(dot.plot.df$Cell.Type)
dot.plot.df$Protein <- as.factor(dot.plot.df$Protein)
dot.plot.df$PES <- abs(dot.plot.df$PES)
## calculate significance
pval.vec <- pnorm(abs(dot.plot.df$NES), lower.tail = FALSE)
pval.vec <- p.adjust(pval.vec, method = 'BH')
names(pval.vec) <- rownames(dot.plot.df)
sig.vec <- rep('NS', length(pval.vec))
names(sig.vec) <- rownames(dot.plot.df)
sig.vec[intersect(names(which(pval.vec < 0.05)), rownames(dot.plot.df)[which(dot.plot.df$NES > 0)])] <- 'UP'
sig.vec[intersect(names(which(pval.vec < 0.05)), rownames(dot.plot.df)[which(dot.plot.df$NES < 0)])] <- 'DOWN'
dot.plot.df$SIG <- as.factor(sig.vec)
## make plot
sig.plot.mat <- dot.plot.df[which(dot.plot.df$SIG != 'NS'),]
ns.plot.mat <- dot.plot.df[which(dot.plot.df$SIG == 'NS'),]
ggplot(sig.plot.mat, aes(x = Protein, y = Cell.Type)) +
  geom_point(aes(size = PES, color = NES)) +
  labs(title = 'WNT TF Activity in 926 vs Vehicle', y = 'Cell Type')   +
  #scale_color_manual(values = sig.cols) + 
  scale_color_gradient2(low = 'blue', mid = 'grey', high = 'red', midpoint = 0) +
  scale_size_continuous(range = c(3, 10), limits = c(0, 0.2)) + 
  guides(size = guide_legend(title="|PES|", order = 2),
         color = guide_colourbar(title = "NES", order = 1)) +
  theme_bw() +
  geom_point(data = ns.plot.mat, aes(Protein, Cell.Type, color = NES, size = abs(PES)), color = 'black')
ggsave('anchor-combined_analysis/anchor_analysis/fig-7_wnt-tf-dif-activity.pdf',
       height = 4, width = 10, units = 'in', dpi = 300)
###############

### 8 - dot plot for endothelial diffeential activity
###############
narnea.obj <- readRDS('anchor-combined_analysis/ct_group-narnea/Endothelial_narnea.rds')
## prepare plot data
use.prots <- intersect(rownames(narnea.obj$NES), endo.genes)
nes.vec <- narnea.obj$NES[use.prots,]
pes.vec <- narnea.obj$PES[use.prots,]
# convert to p-values
pval.vec <- pnorm(abs(nes.vec), lower.tail = FALSE)
sig.vec <- rep('NS', length(use.prots)); names(sig.vec) <- use.prots
sig.vec[intersect(which(pval.vec < 0.05), which(pes.vec > 0))] <- 'UP'
sig.vec[intersect(which(pval.vec < 0.05), which(pes.vec < 0))] <- 'DOWN'
# make plot data
plot.df <- data.frame('sig' = sig.vec,
                      'pes' = abs(pes.vec),
                      'nes' = nes.vec,
                      'gene' = use.prots,
                      'ct' = 'Endothelial')
# make plot
sig.plot.mat <- plot.df[which(plot.df$sig != 'NS'),]
ns.plot.mat <- plot.df[which(plot.df$sig == 'NS'),]
ggplot(sig.plot.mat, aes(gene, ct)) + geom_point(aes(color = nes, size = abs(pes))) +
  labs(x = 'Protein', y = 'Cell Type', title = 'Endothelial Protein Differential Activity') +
  #scale_color_manual(values = sig.cols) + 
  scale_color_gradient2(low = 'blue', mid = 'grey', high = 'red', midpoint = 0) +
  scale_size_continuous(range = c(3, 10), limits = c(0, 0.2)) + 
  guides(size = guide_legend(title="|PES|", order = 2),
         color = guide_colourbar(title = "NES", order = 1))
ggsave('anchor-combined_analysis/anchor_analysis/fig-8_endothelial-dif-activity-dot-plot.pdf',
       height = 4, width = 10, units = 'in', dpi = 300)
###############

### 9 - dot plot for change in WNT expression
###############
seurat.obj <- readRDS('anchor-combined_analysis/anchored_seurat-obj.rds')
cell.names <- colnames(seurat.obj)
cell.samples <- sapply(cell.names, function(x) {strsplit(x, '\\.')[[1]][1]} )
cell.treatment <- sample.table$Group[match(cell.samples, sample.table$Sample)]
names(cell.treatment) <- names(cell.samples)
clust.vec <- seurat.obj$seurat_clusters
ct.vec <- plyr::mapvalues(clust.vec, from = cluster.mappings$Cluster, to = cluster.mappings$Cell.Type)
## load count matrix
cpm.mat <- readRDS('anchor-combined_analysis/combined-gexp_pad-filt-counts_cpm.rds')
## dif gene exp analysis for each cell type
use.wnt.genes <- intersect(wnt.genes, rownames(cpm.mat))
dif.exp.mats <- list()
for (ct in unique(ct.vec)) {
  print(ct)
  # get group labels
  ct.cells <- names(ct.vec)[which(ct.vec == ct)]
  ct.926.cells <- intersect(ct.cells, names(cell.treatment)[which(cell.treatment == '926')])
  ct.vehicle.cells <- intersect(ct.cells, names(cell.treatment)[which(cell.treatment == 'Vehicle')])
  sct.mat <- as.matrix(seurat.obj@assays$SCT@scale.data)
  # perform analysis
  dif.exp <- sapply(use.wnt.genes, function(x) {
    test.vec <- cpm.mat[x, ct.926.cells]
    ref.vec <- cpm.mat[x, ct.vehicle.cells]
    test.obj <- wilcox.test(test.vec, ref.vec, alternative = "two.sided")
    p.val <- test.obj$p.val
    rbs.cor <- 2 * test.obj$statistic / (length(ct.926.cells) * length(ct.vehicle.cells)) - 1
    return(c(x, p.val, rbs.cor))
  })
  # reformat to matrix
  rownames(dif.exp) <- c('gene', 'pval', 'rbsc')
  dif.exp.mat <- as.data.frame(t(dif.exp))
  dif.exp.mat$ct <- ct
  dif.exp.mat$pval <- as.numeric(dif.exp.mat$pval)
  dif.exp.mat$rbsc <- as.numeric(dif.exp.mat$rbsc)
  dif.exp.mat$logp <- log(dif.exp.mat$pval)
  sig.vec <- rep('NS', nrow(dif.exp.mat))
  sig.vec[which(dif.exp.mat$pval < 0.05)] <- 'DOWN'
  dif.exp.mat$sig <- sig.vec
  # calculate NES
  nes.vec <- qnorm(dif.exp.mat$logp, log.p = TRUE, lower.tail = FALSE) * sign(dif.exp.mat$rbsc)
  dif.exp.mat$nes <- nes.vec
  # add to list
  dif.exp.mats[[ct]] <- dif.exp.mat
}
## all cell types
plot.mat <- Reduce(rbind, dif.exp.mats)
sig.plot.mat <- plot.mat[which(plot.mat$sig != 'NS'),]
ns.plot.mat <- plot.mat[which(plot.mat$sig == 'NS'),]
ggplot(sig.plot.mat, aes(gene, ct)) + geom_point(aes(color = nes, size = abs(rbsc))) +
  labs(x = 'Gene', y = 'Cell Type', title = 'WNT Genes Differential Expression - All Cell Types') +
  #scale_color_manual(values = sig.cols) + 
  scale_color_gradient2(low = 'blue', mid = 'grey', high = 'red', midpoint = 0) +
  scale_size_continuous(range = c(3, 10), limits = c(0, 0.2)) + 
  guides(size = guide_legend(title="|RBSC|", order = 2),
         color = guide_colourbar(title = "NES", order = 1)) + theme_bw() +
  geom_point(data = ns.plot.mat, aes(gene, ct, color = nes, size = abs(rbsc)), color = 'black')
ggsave('anchor-combined_analysis/anchor_analysis/fig-9_wnt-genes-all-cells-dif-exp-dot-plot.pdf',
       height = 4, width = 10, units = 'in', dpi = 300)
###############

### 10 - dot plot for change in WNT activity
###############
ct.wnt.dfs <- list()
for (ct in sort(unique(ct.vec))) {
  cat(paste(ct, '\n'))
  narnea.res <- readRDS(paste('anchor-combined_analysis/ct_group-narnea/', ct, '_narnea.rds', sep = ''))
  tf.intersect <- intersect(rownames(narnea.res$NES), wnt.genes)
  tf.nes <- narnea.res$NES[tf.intersect,]
  tf.pes <- narnea.res$PES[tf.intersect,]
  ct.df <- data.frame('Cell.Type' = ct,
                      'Protein' = tf.intersect,
                      'PES' = tf.pes,
                      'NES' = tf.nes)
  ct.wnt.dfs[[ct]] <- ct.df
}
dot.plot.df <- Reduce(rbind, ct.wnt.dfs)
## reformat dataframe
dot.plot.df <- as.data.frame(dot.plot.df)
dot.plot.df$Cell.Type <- as.factor(dot.plot.df$Cell.Type)
dot.plot.df$Protein <- as.factor(dot.plot.df$Protein)
dot.plot.df$PES <- abs(dot.plot.df$PES)
## calculate significance
pval.vec <- pnorm(abs(dot.plot.df$NES), lower.tail = FALSE)
pval.vec <- p.adjust(pval.vec, method = 'BH')
names(pval.vec) <- rownames(dot.plot.df)
sig.vec <- rep('NS', length(pval.vec))
names(sig.vec) <- rownames(dot.plot.df)
sig.vec[intersect(names(which(pval.vec < 0.05)), rownames(dot.plot.df)[which(dot.plot.df$NES > 0)])] <- 'UP'
sig.vec[intersect(names(which(pval.vec < 0.05)), rownames(dot.plot.df)[which(dot.plot.df$NES < 0)])] <- 'DOWN'
dot.plot.df$SIG <- as.factor(sig.vec)
## make plot
sig.plot.mat <- dot.plot.df[which(dot.plot.df$SIG != 'NS'),]
ns.plot.mat <- dot.plot.df[which(dot.plot.df$SIG == 'NS'),]
ggplot(sig.plot.mat, aes(x = Protein, y = Cell.Type)) +
  geom_point(aes(size = PES, color = NES)) +
  labs(title = 'WNT Activity in 926 vs Vehicle', y = 'Cell Type')  +
  #scale_color_manual(values = sig.cols) + 
  scale_color_gradient2(low = 'blue', mid = 'grey', high = 'red', midpoint = 0) +
  scale_size_continuous(range = c(3, 10), limits = c(0, 0.2)) + 
  guides(size = guide_legend(title="|PES|", order = 2),
         color = guide_colourbar(title = "NES", order = 1)) +
  theme_bw() +
  geom_point(data = ns.plot.mat, aes(Protein, Cell.Type, color = NES, size = abs(PES)), color = 'black')
ggsave('anchor-combined_analysis/anchor_analysis/fig-10_wnt-gene-dif-activity.pdf',
       height = 4, width = 10, units = 'in', dpi = 300)
###############

### 11 - dot plot for angio receptor activity
###############
ct.plot.dfs <- list()
for (ct in names(ct.cols)) {
  print(ct)
  narnea.obj <- readRDS(paste('anchor-combined_analysis/ct_group-narnea/', ct, '_narnea.rds', sep = ''))
  # get vectors
  use.prots <- intersect(rownames(narnea.obj$NES), angio.ligands)
  nes.vec <- narnea.obj$NES[use.prots,]
  pes.vec <- narnea.obj$PES[use.prots,]
  # convert to significance
  pval.vec <- pnorm(abs(nes.vec), lower.tail = FALSE)
  sig.vec <- rep('NS', length(use.prots)); names(sig.vec) <- use.prots
  sig.vec[intersect(which(pval.vec < 0.05), which(pes.vec > 0))] <- 'UP'
  sig.vec[intersect(which(pval.vec < 0.05), which(pes.vec < 0))] <- 'DOWN'
  # make plot data
  plot.df <- data.frame('sig' = sig.vec,
                        'pes' = abs(pes.vec),
                        'nes' = nes.vec,
                        'gene' = use.prots,
                        'ct' = ct)
  ct.plot.dfs[[ct]] <- plot.df
}
## make plot
plot.df <- Reduce(rbind, ct.plot.dfs)
sig.plot.mat <- plot.df[which(plot.df$sig != 'NS'),]
ns.plot.mat <- plot.df[which(plot.df$sig == 'NS'),]
ggplot(sig.plot.mat, aes(gene, ct)) + geom_point(aes(color = nes, size = abs(pes))) +
  labs(x = 'Protein', y = 'Cell Type', title = 'Angiogenesis Ligand Protein Differential Activity - All Cell Types') +
  #scale_color_manual(values = sig.cols) + 
  scale_color_gradient2(low = 'blue', mid = 'grey', high = 'red', midpoint = 0) +
  scale_size_continuous(range = c(3, 10), limits = c(0, 0.2)) + 
  guides(size = guide_legend(title="|PES|", order = 2),
         color = guide_colourbar(title = "NES", order = 1)) +
  theme_bw() +
  geom_point(data = ns.plot.mat, aes(gene, ct, color = nes, size = abs(pes)), color = 'black')
ggsave('anchor-combined_analysis/anchor_analysis/fig-11_angio-ligand-dif-activity-dot-plot.pdf',
       height = 4,  width = 10, units = 'in', dpi = 300)
###############

### 12 - dot plot for angio ligand activity
###############
ct.plot.dfs <- list()
for (ct in names(ct.cols)) {
  print(ct)
  narnea.obj <- readRDS(paste('anchor-combined_analysis/ct_group-narnea/', ct, '_narnea.rds', sep = ''))
  # get vectors
  use.prots <- intersect(rownames(narnea.obj$NES), angio.receptors)
  nes.vec <- narnea.obj$NES[use.prots,]
  pes.vec <- narnea.obj$PES[use.prots,]
  # convert to significance
  pval.vec <- pnorm(abs(nes.vec), lower.tail = FALSE)
  sig.vec <- rep('NS', length(use.prots)); names(sig.vec) <- use.prots
  sig.vec[intersect(which(pval.vec < 0.05), which(pes.vec > 0))] <- 'UP'
  sig.vec[intersect(which(pval.vec < 0.05), which(pes.vec < 0))] <- 'DOWN'
  # make plot data
  plot.df <- data.frame('sig' = sig.vec,
                        'pes' = abs(pes.vec),
                        'nes' = nes.vec,
                        'gene' = use.prots,
                        'ct' = ct)
  ct.plot.dfs[[ct]] <- plot.df
}
## make plot
plot.df <- Reduce(rbind, ct.plot.dfs)
sig.plot.mat <- plot.df[which(plot.df$sig != 'NS'),]
ns.plot.mat <- plot.df[which(plot.df$sig == 'NS'),]
ggplot(plot.df, aes(gene, ct)) + geom_point(aes(color = nes, size = abs(pes))) +
  labs(x = 'Protein', y = 'Cell Type', title = 'Angiogenesis Receptors Protein Differential Activity - All Cell Types') +
  #scale_color_manual(values = sig.cols) + 
  scale_color_gradient2(low = 'blue', mid = 'grey', high = 'red', midpoint = 0) +
  scale_size_continuous(range = c(3, 10), limits = c(0, 0.2)) + 
  guides(size = guide_legend(title="|PES|", order = 2),
         color = guide_colourbar(title = "NES", order = 1)) +
  theme_bw() +
  geom_point(data = ns.plot.mat, aes(gene, ct, color = nes, size = abs(pes)), color = 'black')
ggsave('anchor-combined_analysis/anchor_analysis/fig-12_angio-receptors-dif-activity-dot-plot.pdf',
       height = 4,  width = 10, units = 'in', dpi = 300)
###############

### 13 - dot plot for angio receptor expression
###############
seurat.obj <- readRDS('anchor-combined_analysis/anchored_seurat-obj.rds')
cell.names <- colnames(seurat.obj)
cell.samples <- sapply(cell.names, function(x) {strsplit(x, '\\.')[[1]][1]} )
cell.treatment <- sample.table$Group[match(cell.samples, sample.table$Sample)]
names(cell.treatment) <- names(cell.samples)
clust.vec <- seurat.obj$seurat_clusters
ct.vec <- plyr::mapvalues(clust.vec, from = cluster.mappings$Cluster, to = cluster.mappings$Cell.Type)
## load count matrix
cpm.mat <- readRDS('anchor-combined_analysis/combined-gexp_pad-filt-counts_cpm.rds')
## dif gene exp analysis for each cell type
use.genes <- intersect(angio.receptors, rownames(cpm.mat))
dif.exp.mats <- list()
for (ct in unique(ct.vec)) {
  print(ct)
  # get group labels
  ct.cells <- names(ct.vec)[which(ct.vec == ct)]
  ct.926.cells <- intersect(ct.cells, names(cell.treatment)[which(cell.treatment == '926')])
  ct.vehicle.cells <- intersect(ct.cells, names(cell.treatment)[which(cell.treatment == 'Vehicle')])
  sct.mat <- as.matrix(seurat.obj@assays$SCT@scale.data)
  # perform analysis
  dif.exp <- sapply(use.genes, function(x) {
    test.vec <- cpm.mat[x, ct.926.cells]
    ref.vec <- cpm.mat[x, ct.vehicle.cells]
    test.obj <- wilcox.test(test.vec, ref.vec, alternative = "two.sided")
    p.val <- test.obj$p.val
    rbs.cor <- 2 * test.obj$statistic / (length(ct.926.cells) * length(ct.vehicle.cells)) - 1
    return(c(x, p.val, rbs.cor))
  })
  # reformat to matrix
  rownames(dif.exp) <- c('gene', 'pval', 'rbsc')
  dif.exp.mat <- as.data.frame(t(dif.exp))
  dif.exp.mat$ct <- ct
  dif.exp.mat$pval <- as.numeric(dif.exp.mat$pval)
  dif.exp.mat$rbsc <- as.numeric(dif.exp.mat$rbsc)
  dif.exp.mat$logp <- log(dif.exp.mat$pval)
  sig.vec <- rep('NS', nrow(dif.exp.mat))
  sig.vec[which(dif.exp.mat$pval < 0.05)] <- 'DOWN'
  dif.exp.mat$sig <- sig.vec
  # calculate NES
  nes.vec <- qnorm(dif.exp.mat$logp, log.p = TRUE, lower.tail = FALSE) * sign(dif.exp.mat$rbsc)
  dif.exp.mat$nes <- nes.vec
  # add to list
  dif.exp.mats[[ct]] <- dif.exp.mat
}
## all cell types
plot.mat <- Reduce(rbind, dif.exp.mats)
sig.plot.mat <- plot.mat[which(plot.mat$sig != 'NS'),]
ns.plot.mat <- plot.mat[which(plot.mat$sig == 'NS'),]
ggplot(sig.plot.mat, aes(gene, ct)) + geom_point(aes(color = nes, size = abs(rbsc))) +
  labs(x = 'Gene', y = 'Cell Type', title = 'Angiogensis Receptor Genes Differential Expression - All Cell Types') +
  #scale_color_manual(values = sig.cols) + 
  scale_color_gradient2(low = 'blue', mid = 'grey', high = 'red', midpoint = 0) +
  scale_size_continuous(range = c(3, 10), limits = c(0, 0.2)) + 
  guides(size = guide_legend(title="|RBSC|", order = 2),
         color = guide_colourbar(title = "NES", order = 1)) + 
  theme_bw() +
  geom_point(data = ns.plot.mat, aes(gene, ct, color = nes, size = abs(rbsc)), color = 'black')
ggsave('anchor-combined_analysis/anchor_analysis/fig-13_angio-receptors-all-cells-dif-exp-dot-plot.pdf',
       height = 4, width = 10, units = 'in', dpi = 300)
###############

### 14 - dot plot for angio ligand expression
###############
seurat.obj <- readRDS('anchor-combined_analysis/anchored_seurat-obj.rds')
cell.names <- colnames(seurat.obj)
cell.samples <- sapply(cell.names, function(x) {strsplit(x, '\\.')[[1]][1]} )
cell.treatment <- sample.table$Group[match(cell.samples, sample.table$Sample)]
names(cell.treatment) <- names(cell.samples)
clust.vec <- seurat.obj$seurat_clusters
ct.vec <- plyr::mapvalues(clust.vec, from = cluster.mappings$Cluster, to = cluster.mappings$Cell.Type)
## load count matrix
cpm.mat <- readRDS('anchor-combined_analysis/combined-gexp_pad-filt-counts_cpm.rds')
## dif gene exp analysis for each cell type
use.genes <- intersect(angio.ligands, rownames(cpm.mat))
dif.exp.mats <- list()
for (ct in unique(ct.vec)) {
  print(ct)
  # get group labels
  ct.cells <- names(ct.vec)[which(ct.vec == ct)]
  ct.926.cells <- intersect(ct.cells, names(cell.treatment)[which(cell.treatment == '926')])
  ct.vehicle.cells <- intersect(ct.cells, names(cell.treatment)[which(cell.treatment == 'Vehicle')])
  sct.mat <- as.matrix(seurat.obj@assays$SCT@scale.data)
  # perform analysis
  dif.exp <- sapply(use.genes, function(x) {
    test.vec <- cpm.mat[x, ct.926.cells]
    ref.vec <- cpm.mat[x, ct.vehicle.cells]
    test.obj <- wilcox.test(test.vec, ref.vec, alternative = "two.sided")
    p.val <- test.obj$p.val
    rbs.cor <- 2 * test.obj$statistic / (length(ct.926.cells) * length(ct.vehicle.cells)) - 1
    return(c(x, p.val, rbs.cor))
  })
  # reformat to matrix
  rownames(dif.exp) <- c('gene', 'pval', 'rbsc')
  dif.exp.mat <- as.data.frame(t(dif.exp))
  dif.exp.mat$ct <- ct
  dif.exp.mat$pval <- as.numeric(dif.exp.mat$pval)
  dif.exp.mat$rbsc <- as.numeric(dif.exp.mat$rbsc)
  dif.exp.mat$logp <- log(dif.exp.mat$pval)
  sig.vec <- rep('NS', nrow(dif.exp.mat))
  sig.vec[which(dif.exp.mat$pval < 0.05)] <- 'DOWN'
  dif.exp.mat$sig <- sig.vec
  # calculate NES
  nes.vec <- qnorm(dif.exp.mat$logp, log.p = TRUE, lower.tail = FALSE) * sign(dif.exp.mat$rbsc)
  dif.exp.mat$nes <- nes.vec
  # add to list
  dif.exp.mats[[ct]] <- dif.exp.mat
}
## all cell types
plot.mat <- Reduce(rbind, dif.exp.mats)
sig.plot.mat <- plot.mat[which(plot.mat$sig != 'NS'),]
ns.plot.mat <- plot.mat[which(plot.mat$sig == 'NS'),]
ggplot(sig.plot.mat, aes(gene, ct)) + geom_point(aes(color = nes, size = abs(rbsc))) +
  labs(x = 'Gene', y = 'Cell Type', title = 'Angiogensis Ligand Genes Differential Expression - All Cell Types') +
  #scale_color_manual(values = sig.cols) + 
  scale_color_gradient2(low = 'blue', mid = 'grey', high = 'red', midpoint = 0) +
  scale_size_continuous(range = c(3, 10), limits = c(0, 0.2)) + 
  guides(size = guide_legend(title="|RBSC|", order = 2),
         color = guide_colourbar(title = "NES", order = 1)) + 
  theme_bw() +
  geom_point(data = ns.plot.mat, aes(gene, ct, color = nes, size = abs(rbsc)), color = 'black')
ggsave('anchor-combined_analysis/anchor_analysis/fig-14_angio-ligands-all-cells-dif-exp-dot-plot.pdf',
       height = 4, width = 10, units = 'in', dpi = 300)
###############

### 15 - wnt expression in vehicle
###############
seurat.obj <- readRDS('anchor-combined_analysis/anchored_seurat-obj.rds')
cell.names <- colnames(seurat.obj)
cell.samples <- sapply(cell.names, function(x) {strsplit(x, '\\.')[[1]][1]} )
cell.treatment <- sample.table$Group[match(cell.samples, sample.table$Sample)]
names(cell.treatment) <- names(cell.samples)
clust.vec <- seurat.obj$seurat_clusters
ct.vec <- plyr::mapvalues(clust.vec, from = cluster.mappings$Cluster, to = cluster.mappings$Cell.Type)
## load count matrix
cpm.mat <- readRDS('anchor-combined_analysis/combined-gexp_pad-filt-counts_cpm.rds')
vehicle.cells <- names(cell.treatment)[which(cell.treatment == 'Vehicle')]
## get dataframe for each wnt gene
ct.exp.list <- list()
for (gn in intersect(wnt.genes, rownames(cpm.mat))) {
  for (ct in unique(ct.vec)) {
    ct.cells <- intersect(names(ct.vec)[which(ct.vec == ct)], vehicle.cells)
    exp.vec <- cpm.mat[gn, ct.cells]
    mean.exp <- mean(exp.vec)
    exp.per <- length(which(exp.vec > 0)) / length(exp.vec)
    ct.exp.list[[paste(gn, ct, sep = '.')]] <- c(ct, gn, mean.exp, exp.per)
  }
}
wnt.vehicle.exp <- as.data.frame(Reduce(rbind, ct.exp.list))
colnames(wnt.vehicle.exp) <- c('cell.type', 'gene', 'mean.l2cpm', 'exp.per')
wnt.vehicle.exp$mean.l2cpm <- as.numeric(wnt.vehicle.exp$mean.l2cpm)
wnt.vehicle.exp$exp.per <- as.numeric(wnt.vehicle.exp$exp.per)
## make plot
ggplot(wnt.vehicle.exp, aes(gene, cell.type)) + 
  geom_point(aes(color = mean.l2cpm, size = exp.per)) +
  scale_color_gradient(low = 'darkgrey', high = 'red') +
  scale_size_continuous(range = c(3, 10), limits = c(0, 0.2)) + 
  labs(x = 'Gene', y = 'Cell Type', title = 'WNT Ligand Expression in Vehicle Cells') + 
  guides(size = guide_legend(title="% in Group", order = 2),
         color = guide_colourbar(title = "Mean L2CPM", order = 1)) +
  theme_bw()
ggsave('anchor-combined_analysis/anchor_analysis/fig-15_wnt-genes-vehicle-l2cpm.pdf',
       height = 8,  width = 10, units = 'in', dpi = 300)
###############