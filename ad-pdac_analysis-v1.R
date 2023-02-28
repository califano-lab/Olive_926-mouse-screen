setwd('C://Users/Lukas Vlahos/OneDrive/ac_lab/ad_pdac/')
library(PISCES)
library(Seurat)
library(SingleR)
mouse.ref <- MouseRNAseqData()
sample.table <- read.table('sample-sheet.tsv', sep = '\t', header = TRUE)
sample.names <- sample.table$Sample
library(ggplot2)
library(ComplexHeatmap)
library(infercnv)

### read marker csv
###############
ct.markers <- read.csv('pdac_ct-cluster-markers.csv', header = TRUE)
colnames(ct.markers) <- c('Gene', 'CT')
ct.colors <- group_colors(length(unique(ct.markers$CT)), offset = 60)
names(ct.colors) <- unique(ct.markers$CT)
###############

### initial QC
###############
## KA001
s.name <- 'KA001'
# load 10x outs
raw.counts <- Read10X(paste('cellranger_outs/', s.name, '_cellranger_count_outs/filtered_feature_bc_matrix/', sep = ''))
# generate qc plot
qc_plot(as.matrix(raw.counts), species = 'mur', genes = 'symb')
# filter count matrix
filt.counts <- qc_filt(as.matrix(raw.counts), min.depth = 2500, max.depth = 25000, 
                       max.genes = 5000, max.mt = 0.2)
# create directory and save objects
dir.create(paste('sample_analysis/', s.name, sep = ''))
saveRDS(filt.counts, file = paste('sample_analysis/', s.name, '/', s.name, '_filt-counts.rds', sep = ''))
jpeg(filename = paste('sample_analysis/', s.name, '/', s.name, '_qc.jpg', sep = ''),
     height = 6, width = 8, units = 'in', res = 150)
qc_plot(as.matrix(raw.counts), species = 'mur', genes = 'symb')
dev.off()
## KA003
s.name <- 'KA003'
# load 10x outs
raw.counts <- Read10X(paste('cellranger_outs/', s.name, '_cellranger_count_outs/filtered_feature_bc_matrix/', sep = ''))
# generate qc plot
qc_plot(as.matrix(raw.counts), species = 'mur', genes = 'symb')
# filter count matrix
filt.counts <- qc_filt(as.matrix(raw.counts), min.depth = 1400, max.depth = 25000, 
                       max.genes = 6000, max.mt = 0.25)
# create directory and save objects
dir.create(paste('sample_analysis/', s.name, sep = ''))
saveRDS(filt.counts, file = paste('sample_analysis/', s.name, '/', s.name, '_filt-counts.rds', sep = ''))
jpeg(filename = paste('sample_analysis/', s.name, '/', s.name, '_qc.jpg', sep = ''),
     height = 6, width = 8, units = 'in', res = 150)
qc_plot(as.matrix(raw.counts), species = 'mur', genes = 'symb')
dev.off()
## KA005
s.name <- 'KA005'
# load 10x outs
raw.counts <- Read10X(paste('cellranger_outs/', s.name, '_cellranger_count_outs/filtered_feature_bc_matrix/', sep = ''))
# generate qc plot
qc_plot(as.matrix(raw.counts), species = 'mur', genes = 'symb')
# filter count matrix
filt.counts <- qc_filt(as.matrix(raw.counts), min.depth = 1000, max.depth = 25000, 
                       max.genes = 6000, max.mt = 0.25)
# create directory and save objects
dir.create(paste('sample_analysis/', s.name, sep = ''))
saveRDS(filt.counts, file = paste('sample_analysis/', s.name, '/', s.name, '_filt-counts.rds', sep = ''))
jpeg(filename = paste('sample_analysis/', s.name, '/', s.name, '_qc.jpg', sep = ''),
     height = 6, width = 8, units = 'in', res = 150)
qc_plot(as.matrix(raw.counts), species = 'mur', genes = 'symb')
dev.off()
## KA006
s.name <- 'KA006'
# load 10x outs
raw.counts <- Read10X(paste('cellranger_outs/', s.name, '_cellranger_count_outs/filtered_feature_bc_matrix/', sep = ''))
# generate qc plot
qc_plot(as.matrix(raw.counts), species = 'mur', genes = 'symb')
# filter count matrix
filt.counts <- qc_filt(as.matrix(raw.counts), min.depth = 1400, max.depth = 25000, 
                       max.genes = 6000, max.mt = 0.25)
# create directory and save objects
dir.create(paste('sample_analysis/', s.name, sep = ''))
saveRDS(filt.counts, file = paste('sample_analysis/', s.name, '/', s.name, '_filt-counts.rds', sep = ''))
jpeg(filename = paste('sample_analysis/', s.name, '/', s.name, '_qc.jpg', sep = ''),
     height = 6, width = 8, units = 'in', res = 150)
qc_plot(as.matrix(raw.counts), species = 'mur', genes = 'symb')
dev.off()
## KA007
s.name <- 'KA007'
# load 10x outs
raw.counts <- Read10X(paste('cellranger_outs/', s.name, '_cellranger_count_outs/filtered_feature_bc_matrix/', sep = ''))
# generate qc plot
qc_plot(as.matrix(raw.counts), species = 'mur', genes = 'symb')
# filter count matrix
filt.counts <- qc_filt(as.matrix(raw.counts), min.depth = 1500, max.depth = 25000, 
                       max.genes = 6000, max.mt = 0.25)
# create directory and save objects
dir.create(paste('sample_analysis/', s.name, sep = ''))
saveRDS(filt.counts, file = paste('sample_analysis/', s.name, '/', s.name, '_filt-counts.rds', sep = ''))
jpeg(filename = paste('sample_analysis/', s.name, '/', s.name, '_qc.jpg', sep = ''),
     height = 6, width = 8, units = 'in', res = 150)
qc_plot(as.matrix(raw.counts), species = 'mur', genes = 'symb')
dev.off()
###############

### process Lorenzo vehicle samples
###############
saline.sur <- 'sal'
## KA001
sn <- 'KA001'
seurat.obj <- readRDS('prev_vehicle/KA001_Seurat.rds')
counts.mat <- as.matrix(seurat.obj@assays$RNA@counts)
qc_plot(counts.mat, species = 'mur', genes = 'symb')
# create directory and save objects
s.name <- paste(sn, saline.sur, sep = '-')
dir.create(paste('sample_analysis/', s.name, sep = ''))
saveRDS(counts.mat, file = paste('sample_analysis/', s.name, '/', s.name, '_filt-counts.rds', sep = ''))
jpeg(filename = paste('sample_analysis/', s.name, '/', s.name, '_qc.jpg', sep = ''),
     height = 6, width = 8, units = 'in', res = 150)
qc_plot(as.matrix(counts.mat), species = 'mur', genes = 'symb')
dev.off()
## KA004
sn <- 'KA004'
seurat.obj <- readRDS('prev_vehicle/KA001_Seurat.rds')
counts.mat <- as.matrix(seurat.obj@assays$RNA@counts)
qc_plot(counts.mat, species = 'mur', genes = 'symb')
# create directory and save objects
s.name <- paste(sn, saline.sur, sep = '-')
dir.create(paste('sample_analysis/', s.name, sep = ''))
saveRDS(counts.mat, file = paste('sample_analysis/', s.name, '/', s.name, '_filt-counts.rds', sep = ''))
jpeg(filename = paste('sample_analysis/', s.name, '/', s.name, '_qc.jpg', sep = ''),
     height = 6, width = 8, units = 'in', res = 150)
qc_plot(as.matrix(counts.mat), species = 'mur', genes = 'symb')
dev.off()
## KU030
sn <- 'KA030'
seurat.obj <- readRDS('prev_vehicle/KA001_Seurat.rds')
counts.mat <- as.matrix(seurat.obj@assays$RNA@counts)
qc_plot(counts.mat, species = 'mur', genes = 'symb')
# create directory and save objects
s.name <- paste(sn, saline.sur, sep = '-')
dir.create(paste('sample_analysis/', s.name, sep = ''))
saveRDS(counts.mat, file = paste('sample_analysis/', s.name, '/', s.name, '_filt-counts.rds', sep = ''))
jpeg(filename = paste('sample_analysis/', s.name, '/', s.name, '_qc.jpg', sep = ''),
     height = 6, width = 8, units = 'in', res = 150)
qc_plot(as.matrix(counts.mat), species = 'mur', genes = 'symb')
dev.off()
## KU034
sn <- 'KA034'
seurat.obj <- readRDS('prev_vehicle/KA001_Seurat.rds')
counts.mat <- as.matrix(seurat.obj@assays$RNA@counts)
qc_plot(counts.mat, species = 'mur', genes = 'symb')
# create directory and save objects
s.name <- paste(sn, saline.sur, sep = '-')
dir.create(paste('sample_analysis/', s.name, sep = ''))
saveRDS(counts.mat, file = paste('sample_analysis/', s.name, '/', s.name, '_filt-counts.rds', sep = ''))
jpeg(filename = paste('sample_analysis/', s.name, '/', s.name, '_qc.jpg', sep = ''),
     height = 6, width = 8, units = 'in', res = 150)
qc_plot(as.matrix(counts.mat), species = 'mur', genes = 'symb')
dev.off()
###############

### gexp PCA plot
###############
## load count matrices and combine
sample.counts <- list()
for (sn in sample.names) {
  c.mat <- readRDS(paste('sample_analysis/', sn, '/', sn, '_filt-counts.rds', sep = ''))
  sample.counts[[sn]] <- c.mat
}
shared.genes <- Reduce(intersect, lapply(sample.counts, function(x) {rownames(x)}))
counts.mat <- Reduce(cbind, lapply(sample.counts, function(x) {x[shared.genes,]}))
sample.name.vec <- rep(sample.names, sapply(sample.counts, ncol))
sample.group.vec <- rep(sample.table$Group, sapply(sample.counts, ncol))
## cpm normalize and run pca
norm.mat <- cpm_norm(counts.mat, l2 = TRUE)
pca.obj <- fast_pca(norm.mat)
saveRDS(pca.obj, file = 'combined_analysis/combined-expression_cpm-pca.rds')
## make plot
plot.df <- data.frame('PC1' = pca.obj$x[,1], 'PC2' = pca.obj$x[,2],
                      'Sample' = sample.name.vec,
                      'Group' = sample.group.vec)
# sample
ggplot(plot.df, aes(PC1, PC2)) + geom_point(aes(color = Sample), size = 0.5) +
  facet_wrap(vars(Sample)) +
  ggtitle('Combined Expression PCA')
ggsave('combined_analysis/combined-expression_pca_sample.jpg', 
       height = 6, width = 8, dpi = 300, units = 'in')
# group
ggplot(plot.df, aes(PC1, PC2)) + geom_point(aes(color = Group), size = 0.5) +
  facet_wrap(vars(Group), ncol = 2) +
  ggtitle('Combined Expression PCA')
ggsave('combined_analysis/combined-expression_pca_group.jpg', 
       height = 6, width = 8, dpi = 300, units = 'in')
###############

### gexp clustering w/out saline
###############
sample.inds <- which(sample.table$Group != 'Saline')
sample.names <- sample.table$Sample[sample.inds]
## load count matrices and combine
sample.counts <- list()
for (sn in sample.names) {
  print(sn)
  c.mat <- readRDS(paste('sample_analysis/', sn, '/', sn, '_filt-counts.rds', sep = ''))
  sample.counts[[sn]] <- c.mat
}
shared.genes <- Reduce(intersect, lapply(sample.counts, function(x) {rownames(x)}))
counts.mat <- Reduce(cbind, lapply(sample.counts, function(x) {x[shared.genes,]}))
sample.name.vec <- rep(sample.names, sapply(sample.counts, ncol))
sample.group.vec <- rep(sample.table$Group[sample.inds], sapply(sample.counts, ncol))
colnames(counts.mat) <- paste(sample.name.vec, colnames(counts.mat), sep = '.')
saveRDS(counts.mat, file = 'combined_analysis/combined-gexp_filt-counts.rds')
## create seurat object and cluster
seurat.obj <- CreateSeuratObject(counts = counts.mat, assay = 'RNA', project = 'combo.seur')
seurat.obj <- PercentageFeatureSet(seurat.obj, pattern = "^mt-", col.name = "percent.mt")
seurat.obj <- SCTransform(seurat.obj, vars.to.regress = "percent.mt", verbose = FALSE)
seurat.obj <- RunPCA(seurat.obj, verbose = FALSE)
seurat.obj <- RunUMAP(seurat.obj, dims = 1:30, verbose = FALSE)
seurat.obj <- FindNeighbors(seurat.obj, dims = 1:30, verbose = FALSE)
seurat.obj <- FindClusters(seurat.obj, verbose = FALSE)
saveRDS(seurat.obj, file = 'combined_analysis/combined-gexp_seurat-obj.rds')
## seurat clust markers
cluster.markers <- list()
for (clust.name in levels(seurat.obj$seurat_clusters)) {
  print(clust.name)
  cn.markers <- FindMarkers(seurat.obj, ident.1 = clust.name)
  cn.markers <- cn.markers[order(cn.markers$p_val_adj, -abs(cn.markers$avg_log2FC)),]
  cluster.markers[[as.character(clust.name)]] <- cn.markers
}
saveRDS(cluster.markers, file = 'combined_analysis/combined-gexp_seurat-clust-markers.rds')
# save table of top and bottom markers
top.markers <- lapply(cluster.markers, function(x) {
  c.markers <- x[order(x$avg_log2FC, decreasing = TRUE),]
  return(data.frame('gene' = rownames(c.markers)[1:50],
                    'log2FC' = c.markers$avg_log2FC[1:50]))
})
top.marker.df <- Reduce(cbind, top.markers)
colnames(top.marker.df) <- paste(rep(paste('c', names(top.markers), sep = ''), each = 2),
                                 colnames(top.marker.df), sep = '.')
write.table(top.marker.df, file = 'combined_analysis/combined-gexp_seurat-clust_top-markers.csv',
            sep = ',', col.names = TRUE, row.names = FALSE, quote = FALSE)
bot.markers <- lapply(cluster.markers, function(x) {
  c.markers <- x[order(x$avg_log2FC, decreasing = FALSE),]
  return(data.frame('gene' = rownames(c.markers)[1:50],
                    'log2FC' = c.markers$avg_log2FC[1:50]))
})
bot.marker.df <- Reduce(cbind, bot.markers)
colnames(bot.marker.df) <- paste(rep(paste('c', names(bot.markers), sep = ''), each = 2),
                                 colnames(bot.marker.df), sep = '.')
write.table(bot.marker.df, file = 'combined_analysis/combined-gexp_seurat-clust_bot-markers.csv',
            sep = ',', col.names = TRUE, row.names = FALSE, quote = FALSE)
## singleR
singleR.obj <- SingleR(test = counts.mat, ref = mouse.ref,
                       labels = mouse.ref$label.main, assay.type.test = 1)
saveRDS(singleR.obj, file = 'combined_analysis/combined-gexp_singleR.rds')
singleR.vec <- singleR.obj$first.labels
names(singleR.vec) <- colnames(counts.mat)
###############

### gexp clustering v3
###############
sample.names <- c('KA003', 'KA005', 'KA006', 'KA007')
sample.inds <- match(sample.names, sample.table$Sample)
## load count matrices and combine
sample.counts <- list()
for (sn in sample.names) {
  print(sn)
  c.mat <- readRDS(paste('sample_analysis/', sn, '/', sn, '_filt-counts.rds', sep = ''))
  sample.counts[[sn]] <- c.mat
}
shared.genes <- Reduce(intersect, lapply(sample.counts, function(x) {rownames(x)}))
counts.mat <- Reduce(cbind, lapply(sample.counts, function(x) {x[shared.genes,]}))
sample.name.vec <- rep(sample.names, sapply(sample.counts, ncol))
sample.group.vec <- rep(sample.table$Group[sample.inds], sapply(sample.counts, ncol))
names(sample.name.vec) <- colnames(counts.mat)
names(sample.group.vec) <- colnames(counts.mat)
colnames(counts.mat) <- paste(sample.name.vec, colnames(counts.mat), sep = '.')
saveRDS(counts.mat, file = 'combined_analysis/combined-gexp_filt-counts.rds')
## create seurat object and cluster
seurat.obj <- CreateSeuratObject(counts = counts.mat, assay = 'RNA', project = 'combo.seur')
seurat.obj <- PercentageFeatureSet(seurat.obj, pattern = "^mt-", col.name = "percent.mt")
seurat.obj <- SCTransform(seurat.obj, vars.to.regress = "percent.mt", verbose = FALSE)
seurat.obj <- RunPCA(seurat.obj, verbose = FALSE)
seurat.obj <- RunUMAP(seurat.obj, dims = 1:30, verbose = FALSE)
seurat.obj <- FindNeighbors(seurat.obj, dims = 1:30, verbose = FALSE)
seurat.obj <- FindClusters(seurat.obj, verbose = FALSE)
saveRDS(seurat.obj, file = 'combined_analysis/combined-gexp_seurat-obj.rds')
## seurat clust markers
cluster.markers <- list()
for (clust.name in levels(seurat.obj$seurat_clusters)) {
  print(clust.name)
  cn.markers <- FindMarkers(seurat.obj, ident.1 = clust.name)
  cn.markers <- cn.markers[order(cn.markers$p_val_adj, -abs(cn.markers$avg_log2FC)),]
  cluster.markers[[as.character(clust.name)]] <- cn.markers
}
saveRDS(cluster.markers, file = 'combined_analysis/combined-gexp_seurat-clust-markers.rds')
# save table of top and bottom markers
top.markers <- lapply(cluster.markers, function(x) {
  c.markers <- x[order(x$avg_log2FC, decreasing = TRUE),]
  return(data.frame('gene' = rownames(c.markers)[1:50],
                    'log2FC' = c.markers$avg_log2FC[1:50]))
})
top.marker.df <- Reduce(cbind, top.markers)
colnames(top.marker.df) <- paste(rep(paste('c', names(top.markers), sep = ''), each = 2),
                                 colnames(top.marker.df), sep = '.')
write.table(top.marker.df, file = 'combined_analysis/combined-gexp_seurat-clust_top-markers.csv',
            sep = ',', col.names = TRUE, row.names = FALSE, quote = FALSE)
bot.markers <- lapply(cluster.markers, function(x) {
  c.markers <- x[order(x$avg_log2FC, decreasing = FALSE),]
  return(data.frame('gene' = rownames(c.markers)[1:50],
                    'log2FC' = c.markers$avg_log2FC[1:50]))
})
bot.marker.df <- Reduce(cbind, bot.markers)
colnames(bot.marker.df) <- paste(rep(paste('c', names(bot.markers), sep = ''), each = 2),
                                 colnames(bot.marker.df), sep = '.')
write.table(bot.marker.df, file = 'combined_analysis/combined-gexp_seurat-clust_bot-markers.csv',
            sep = ',', col.names = TRUE, row.names = FALSE, quote = FALSE)
## singleR
singleR.obj <- SingleR(test = counts.mat, ref = mouse.ref,
                       labels = mouse.ref$label.main, assay.type.test = 1)
saveRDS(singleR.obj, file = 'combined_analysis/combined-gexp_singleR.rds')
singleR.vec <- singleR.obj$first.labels
names(singleR.vec) <- colnames(counts.mat)
## umap plots
plot.df <- data.frame('UMAP1' = seurat.obj@reductions$umap@cell.embeddings[,1],
                      'UMAP2' = seurat.obj@reductions$umap@cell.embeddings[,2],
                      'Sample' = sample.name.vec, 'Group' = sample.group.vec,
                      'Cluster' = seurat.obj$seurat_clusters,
                      'SingleR' = singleR.vec)
cluster.umap.centers <- lapply(sort(unique(seurat.obj$seurat_clusters)), function(x) {
  clust.samps <- which(seurat.obj$seurat_clusters == x)
  mean.vals <- colMeans(seurat.obj@reductions$umap@cell.embeddings[clust.samps, 1:2])
  return(mean.vals)
})
cluster.umap.centers <- Reduce(rbind, cluster.umap.centers)
cluster.umap.centers <- as.data.frame(cluster.umap.centers)
cluster.umap.centers$Cluster <- sort(unique(seurat.obj$seurat_clusters))
# sample plot
ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = Sample), size = 0.5) +
  ggtitle('Combined Analysis - Sample') +
  guides(colour = guide_legend(override.aes = list(size = 3)))
ggsave('combined_analysis/cluster_plots/combined-gexp_sample-umap.jpg',
       height = 8, width = 10, units = 'in', dpi = 300)
ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = Group), size = 0.5) +
  ggtitle('Combined Analysis - Group') +
  guides(colour = guide_legend(override.aes = list(size = 3)))
ggsave('combined_analysis/cluster_plots/combined-gexp_group-umap.jpg',
       height = 8, width = 10, units = 'in', dpi = 300)
ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = Cluster), size = 0.5) +
  ggtitle('Combined Analysis - Cluster') +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  annotate("text", x = cluster.umap.centers$UMAP_1, y = cluster.umap.centers$UMAP_2,
           label = cluster.umap.centers$Cluster)
ggsave('combined_analysis/cluster_plots/combined-gexp_cluster-umap.jpg',
       height = 8, width = 10, units = 'in', dpi = 300)
ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = SingleR), size = 0.5) +
  ggtitle('Combined Analysis - SingleR') +
  guides(colour = guide_legend(override.aes = list(size = 3)))
ggsave('combined_analysis/cluster_plots/combined-gexp_singler-umap.jpg',
       height = 8, width = 10, units = 'in', dpi = 300)
## marker heatmap
sct.mat <- seurat.obj@assays$SCT@scale.data
cell.order <- names(sort(seurat.obj$seurat_clusters))
marker.use.inds <- which(ct.markers$Gene %in% rownames(sct.mat))
plot.mat <- as.matrix(sct.mat[ct.markers$Gene[marker.use.inds], cell.order])
plot.mat <- t(apply(plot.mat, 1, scale))
# set colors
col.breaks <- quantile_breaks(plot.mat)
col.fun <- circlize::colorRamp2(col.breaks, color_levels('gexp', length(col.breaks)))
clust.colors <- group_colors(length(unique(seurat.obj$seurat_clusters)))
names(clust.colors) <- sort(unique(seurat.obj$seurat_clusters))
singleR.colors <- group_colors(length(unique(singleR.vec)))
names(singleR.colors) <- sort(unique(singleR.vec))
sample.colors <- group_colors(length(unique(sample.name.vec)))
names(sample.colors) <- sort(unique(sample.name.vec))
group.colors <- c('darkgrey', 'lightgrey')
names(group.colors) <- sort(unique(sample.group.vec))
# column annotation
column.annot <- HeatmapAnnotation('Cluster' = seurat.obj$seurat_clusters[cell.order],
                                  'SingleR' = singleR.vec[cell.order],
                                  'Sample' = sample.name.vec[cell.order],
                                  'Group' = sample.group.vec[cell.order],
                                  col = list('Cluster' = clust.colors, 'SingleR' = singleR.colors,
                                             'Sample' = sample.colors, 'Group' = group.colors))
col.gaps <- seurat.obj$seurat_clusters[cell.order]
# row annotation
row.annot <- rowAnnotation('Cell Type' = ct.markers$CT[marker.use.inds],
                           col = list('Cell Type' = ct.colors), 
                           show_annotation_name = FALSE,
                           show_legend = TRUE)
row.gaps <- ct.markers$CT[marker.use.inds]
# make plot
jpeg('combined_analysis/cluster_plots/combined-gexp_cluster-markers.jpg',
     height = 10, width = 15, units = 'in', res = 300)
print(Heatmap(plot.mat, name = 'GEXP',
              col = col.fun,
              top_annotation = column.annot, column_split = col.gaps,
              left_annotation = row.annot, row_split = row.gaps,
              cluster_rows = FALSE, cluster_columns = FALSE,
              show_row_names = TRUE, show_column_names = FALSE,
              column_title = 'Combined Analysis - Supervised Markers', row_title = NULL,
              row_names_gp = gpar(fontsize = 8)))
dev.off()
###############

### full inferCNV analysis
###############
sample.names <- c('KA003', 'KA005', 'KA006', 'KA007')
filt.counts <- readRDS('combined_analysis/combined-gexp_filt-counts.rds')
seurat.obj <- readRDS('combined_analysis/combined-gexp_seurat-obj.rds')
cluster.mappings <- read.csv('combined_analysis/cluster-mappings_20221027.csv')
colnames(cluster.mappings) <- c('Cluster', 'Cell.Type')
## create labeling vector
clust.vec <- seurat.obj$seurat_clusters
ref.clusters <- cluster.mappings$Cluster[which(cluster.mappings$Cell.Type == 'CAFs')]
## write objects as tables
icnv.mat <- filt.counts[, names(clust.vec)]
icnv.ref <- as.character(clust.vec)
names(icnv.ref) <- names(clust.vec)
write.table(icnv.mat, file = 'combined_analysis/test_infercnv/combined_icnv-in.txt', sep = '\t', quote = FALSE)
write.table(icnv.ref, file = 'combined_analysis/test_infercnv/combined_icnv-annot.txt', sep = '\t', quote = FALSE, col.names = FALSE)
## run infercnv
icnv.obj <- CreateInfercnvObject(raw_counts_matrix = 'combined_analysis/test_infercnv/combined_icnv-in.txt',
                                 annotations_file = 'combined_analysis/test_infercnv/combined_icnv-annot.txt',
                                 gene_order_file = '../reference_genomes/mouse_chr_gene_locations.txt',
                                 delim = '\t', ref_group_names = as.character(ref.clusters))
icnv.cells <- run(icnv.obj, cutoff = 0.1, out_dir = 'combined_analysis/test_infercnv/outs/',
                  cluster_by_groups = TRUE, k_obs_groups = 1, analysis_mode = 'samples', 
                  denoise = TRUE, noise_logistic = TRUE,
                  scale_data = FALSE,
                  HMM = FALSE, HMM_type = 'i6', 
                  no_plot = FALSE, plot_steps = TRUE, debug = FALSE, num_threads = 11, up_to_step = 50)
###############

### TEST inferCNV analysis
###############
sample.names <- c('KA003', 'KA005', 'KA006', 'KA007')
filt.counts <- readRDS('combined_analysis/combined-gexp_filt-counts.rds')
seurat.obj <- readRDS('combined_analysis/combined-gexp_seurat-obj.rds')
cluster.mappings <- read.csv('combined_analysis/cluster-mappings_20221027.csv')
colnames(cluster.mappings) <- c('Cluster', 'Cell.Type')
## create labeling vector
clust.vec <- seurat.obj$seurat_clusters
ref.clusters <- cluster.mappings$Cluster[which(cluster.mappings$Cell.Type == 'CAFs')]
## subset to 10 cells from each cluster
subset.cells <- lapply(sort(unique(clust.vec)), function(x) {
  clust.cells <- names(clust.vec)[which(clust.vec == x)]
  clust.sub.cells <- sample(clust.cells, 10)
  return(clust.sub.cells)
})
subset.cells <- unlist(subset.cells)
## write objects as tables
icnv.mat <- filt.counts[, subset.cells]
icnv.ref <- as.character(clust.vec[subset.cells])
names(icnv.ref) <- subset.cells
write.table(icnv.mat, file = 'combined_analysis/test_infercnv/combined_icnv-in.txt', sep = '\t', quote = FALSE)
write.table(icnv.ref, file = 'combined_analysis/test_infercnv/combined_icnv-annot.txt', sep = '\t', quote = FALSE, col.names = FALSE)
## run infercnv
icnv.obj <- CreateInfercnvObject(raw_counts_matrix = 'combined_analysis/test_infercnv/combined_icnv-in.txt',
                                 annotations_file = 'combined_analysis/test_infercnv/combined_icnv-annot.txt',
                                 gene_order_file = '../reference_genomes/mouse_chr_gene_locations.txt',
                                 delim = '\t', ref_group_names = as.character(ref.clusters))
icnv.cells <- run(icnv.obj, cutoff = 0.1, out_dir = 'combined_analysis/test_infercnv/outs/',
                  cluster_by_groups = TRUE, k_obs_groups = 1, analysis_mode = 'samples', 
                  denoise = TRUE, noise_logistic = TRUE,
                  scale_data = FALSE,
                  HMM = FALSE, HMM_type = 'i6', 
                  no_plot = FALSE, plot_steps = TRUE, debug = FALSE, num_threads = 11, up_to_step = 50)
###############

### sample specific inferCNV based on combined analysis clustering
## fibroblast reference
###############
sample.names <- c('KA003', 'KA005', 'KA006', 'KA007')
filt.counts <- readRDS('combined_analysis/combined-gexp_filt-counts.rds')
seurat.obj <- readRDS('combined_analysis/combined-gexp_seurat-obj.rds')
cluster.mappings <- read.csv('combined_analysis/cluster-mappings_20221027.csv')
colnames(cluster.mappings) <- c('Cluster', 'Cell.Type')
sample.name.vec <- sapply(colnames(seurat.obj), function(x) {
  strsplit(x, '\\.')[[1]][1]
})
## create labeling vector
clust.vec <- seurat.obj$seurat_clusters
ref.clusters <- cluster.mappings$Cluster[which(cluster.mappings$Cell.Type == 'CAFs')]
ref.clusters <- as.character(ref.clusters)
test.clusters <- cluster.mappings$Cluster[which(cluster.mappings$Cell.Type %in% c('Epithelial', 'MalEpithelial', 'Undefined'))]
test.clusters <- as.character(test.clusters)
## perform analysis for each sample
for (sn in sample.names) {
  # create directory
  out.dir <- paste('combined_analysis/sample_infercnv/', sn, '_fibro-ref/', sep = '')
  dir.create(out.dir)
  # separate counts and vector
    sn.cells <- names(which(sample.name.vec == sn))
  sn.counts <- filt.counts[, sn.cells]
  sn.clusts <- as.character(clust.vec[sn.cells])
  names(sn.clusts) <- sn.cells
  # subset to the reference and test clusters
  ref.cells <- names(sn.clusts)[which(sn.clusts %in% ref.clusters)]
  test.cells <- names(sn.clusts)[which(sn.clusts %in% test.clusters)]
  icnv.mat <- sn.counts[, c(ref.cells, test.cells)]
  icnv.ref <- sn.clusts[c(ref.cells, test.cells)]
  # save objects
  write.table(icnv.mat, sep = '\t', quote = FALSE,
              file = paste(out.dir, 'combined_icnv-in.txt', sep = ''))
  write.table(icnv.ref, sep = '\t', quote = FALSE, col.names = FALSE,
              file = paste(out.dir, 'combined_icnv-annot.txt', sep = ''))
  # run infercnv
  icnv.obj <- CreateInfercnvObject(raw_counts_matrix = paste(out.dir, 'combined_icnv-in.txt', sep = ''),
                                   annotations_file = paste(out.dir, 'combined_icnv-annot.txt', sep = ''),
                                   gene_order_file = '../reference_genomes/mouse_chr_gene_locations.txt',
                                   delim = '\t', ref_group_names = ref.clusters)
  icnv.cells <- run(icnv.obj, cutoff = 0.1, out_dir = out.dir,
                    cluster_by_groups = TRUE, k_obs_groups = 1, analysis_mode = 'samples', 
                    denoise = TRUE, noise_logistic = TRUE,
                    scale_data = FALSE,
                    HMM = FALSE, HMM_type = 'i6', 
                    no_plot = FALSE, plot_steps = TRUE, debug = FALSE, num_threads = 11, up_to_step = 50)
}
###############

### sample specific inferCNV based on combined analysis clustering
## fibroblast + endothelial reference
###############
sample.names <- c('KA003', 'KA005', 'KA006', 'KA007')
filt.counts <- readRDS('combined_analysis/combined-gexp_filt-counts.rds')
seurat.obj <- readRDS('combined_analysis/combined-gexp_seurat-obj.rds')
cluster.mappings <- read.csv('combined_analysis/cluster-mappings_20221027.csv')
colnames(cluster.mappings) <- c('Cluster', 'Cell.Type')
sample.name.vec <- sapply(colnames(seurat.obj), function(x) {
  strsplit(x, '\\.')[[1]][1]
})
## create labeling vector
clust.vec <- seurat.obj$seurat_clusters
ref.clusters <- cluster.mappings$Cluster[which(cluster.mappings$Cell.Type %in% c('CAFs', 'Endothelial'))]
ref.clusters <- as.character(ref.clusters)
test.clusters <- cluster.mappings$Cluster[which(cluster.mappings$Cell.Type %in% c('Epithelial', 'MalEpithelial', 'Undefined'))]
test.clusters <- as.character(test.clusters)
## perform analysis for each sample
for (sn in sample.names) {
  # create directory
  out.dir <- paste('combined_analysis/sample_infercnv/', sn, '_fibro-endo-ref/', sep = '')
  dir.create(out.dir)
  # separate counts and vector
  sn.cells <- names(which(sample.name.vec == sn))
  sn.counts <- filt.counts[, sn.cells]
  sn.clusts <- as.character(clust.vec[sn.cells])
  names(sn.clusts) <- sn.cells
  # subset to the reference and test clusters
  ref.cells <- names(sn.clusts)[which(sn.clusts %in% ref.clusters)]
  test.cells <- names(sn.clusts)[which(sn.clusts %in% test.clusters)]
  icnv.mat <- sn.counts[, c(ref.cells, test.cells)]
  icnv.ref <- sn.clusts[c(ref.cells, test.cells)]
  # save objects
  write.table(icnv.mat, sep = '\t', quote = FALSE,
              file = paste(out.dir, 'combined_icnv-in.txt', sep = ''))
  write.table(icnv.ref, sep = '\t', quote = FALSE, col.names = FALSE,
              file = paste(out.dir, 'combined_icnv-annot.txt', sep = ''))
  # run infercnv
  icnv.obj <- CreateInfercnvObject(raw_counts_matrix = paste(out.dir, 'combined_icnv-in.txt', sep = ''),
                                   annotations_file = paste(out.dir, 'combined_icnv-annot.txt', sep = ''),
                                   gene_order_file = '../reference_genomes/mouse_chr_gene_locations.txt',
                                   delim = '\t', ref_group_names = ref.clusters)
  icnv.cells <- run(icnv.obj, cutoff = 0.1, out_dir = out.dir,
                    cluster_by_groups = TRUE, k_obs_groups = 1, analysis_mode = 'samples', 
                    denoise = TRUE, noise_logistic = TRUE,
                    scale_data = FALSE,
                    HMM = FALSE, HMM_type = 'i6', 
                    no_plot = FALSE, plot_steps = TRUE, debug = FALSE, num_threads = 11, up_to_step = 50)
}
############### 

### make gene expression plots
###############
seurat.obj <- readRDS('combined_analysis/combined-gexp_seurat-obj.rds')
cluster.markers <- readRDS('combined_analysis/combined-gexp_seurat-clust-markers.rds')
singleR.obj <- readRDS('combined_analysis/combined-gexp_singleR.rds')
singleR.vec <- singleR.obj$first.labels
names(singleR.vec) <- colnames(seurat.obj)
sample.name.vec <- sapply(colnames(seurat.obj), function(x) {
  strsplit(x, '\\.')[[1]][1]
})
sample.group.vec <- sample.table$Group[match(sample.name.vec, sample.table$Sample)]
names(sample.group.vec) <- names(sample.name.vec)
## make scatter plots
plot.df <- data.frame('UMAP1' = seurat.obj@reductions$umap@cell.embeddings[,1],
                      'UMAP2' = seurat.obj@reductions$umap@cell.embeddings[,2],
                      'PCA1' = seurat.obj@reductions$pca@cell.embeddings[,1],
                      'PCA2' = seurat.obj@reductions$pca@cell.embeddings[,2],
                      'Cluster' = seurat.obj$seurat_clusters,
                      'SingleR' = singleR.vec,
                      'Sample' = sample.name.vec,
                      'Group' = sample.group.vec)
# umap plots
ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = Cluster), size = 0.5) +
  ggtitle('Cluster') +
  guides(colour = guide_legend(override.aes = list(size = 3)))
ggsave('combined_analysis/cluster_plots/seurat-clust_umap.jpg',
       height = 8, width = 10, units = 'in', dpi = 300)
ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = SingleR), size = 0.5) +
  ggtitle('SingleR') +
  guides(colour = guide_legend(override.aes = list(size = 3)))
ggsave('combined_analysis/cluster_plots/singleR_umap.jpg',
       height = 8, width = 10, units = 'in', dpi = 300)
ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = Sample), size = 0.5) +
  ggtitle('Sample') +
  guides(colour = guide_legend(override.aes = list(size = 3)))
ggsave('combined_analysis/cluster_plots/sample_umap.jpg',
       height = 8, width = 10, units = 'in', dpi = 300)
ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = Group), size = 0.5) +
  ggtitle('Group') +
  guides(colour = guide_legend(override.aes = list(size = 3)))
ggsave('combined_analysis/cluster_plots/group_umap.jpg',
       height = 8, width = 10, units = 'in', dpi = 300)
# pca plots
ggplot(plot.df, aes(PCA1, PCA2)) + geom_point(aes(color = Sample), size = 0.5) +
  ggtitle('Sample') +
  guides(colour = guide_legend(override.aes = list(size = 3)))
ggsave('combined_analysis/cluster_plots/sample_pca.jpg',
       height = 8, width = 10, units = 'in', dpi = 300)
ggplot(plot.df, aes(PCA1, PCA2)) + geom_point(aes(color = Group), size = 0.5) +
  ggtitle('Group') +
  guides(colour = guide_legend(override.aes = list(size = 3)))
ggsave('combined_analysis/cluster_plots/group_pca.jpg',
       height = 8, width = 10, units = 'in', dpi = 300)
## make heatmap
cell.order <- names(sort(seurat.obj$seurat_clusters))
num.markers <- 10
marker.set <- lapply(names(cluster.markers), function(x) {
  cm.markers <- cluster.markers[[x]]
  c.markers <- cm.markers[order(cm.markers$avg_log2FC, decreasing = TRUE),]
  return(data.frame('gene' = rownames(c.markers)[1:num.markers],
                    'group' = rep(x, num.markers)))
})
marker.set <- Reduce(rbind, marker.set)
marker.set <- marker.set[which(marker.set$gene %in% rownames(seurat.obj@assays$SCT@scale.data)),]
marker.set$group <- as.numeric(marker.set$group)
# set colors
clust.colors <- group_colors(length(unique(seurat.obj$seurat_clusters)))
names(clust.colors) <- sort(unique(seurat.obj$seurat_clusters))
singleR.colors <- group_colors(length(unique(singleR.vec)))
names(singleR.colors) <- sort(unique(singleR.vec))
# create annotations
column.annot <- HeatmapAnnotation('Cluster' = seurat.obj$seurat_clusters[cell.order],
                                  'SingleR' = singleR.vec[cell.order],
                                  'Sample' = sample.name.vec[cell.order],
                                  'Group' = sample.group.vec[cell.order],
                                  col = list('Cluster' = clust.colors, 'SingleR' = singleR.colors))
marker.annot <- rowAnnotation('Cluster' = marker.set[,2],
                              col = list('Cluster' = clust.colors), 
                              show_annotation_name = FALSE,
                              show_legend = FALSE) 
# set gaps
col.gaps <- seurat.obj$seurat_clusters[cell.order]
marker.gaps <- marker.set[,2]
# plot markers
plot.mat <- as.matrix(seurat.obj@assays$SCT@scale.data)[marker.set[,1], cell.order]
plot.mat <- t(apply(plot.mat, 1, scale))
col.breaks <- quantile_breaks(plot.mat)
col.fun <- circlize::colorRamp2(col.breaks, color_levels('gexp', length(col.breaks)))
jpeg('combined_analysis/cluster_plots/cluster-markers_heatmap.jpg',
     height = 20, width = 15, units = 'in', res = 300)
Heatmap(plot.mat, name = 'GEXP',
        col = col.fun,
        top_annotation = column.annot, column_split = col.gaps,
        left_annotation = marker.annot, row_split = marker.gaps,
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = FALSE,
        column_title = 'Seurat Cluster Markers', row_title = NULL,
        row_names_gp = gpar(fontsize = 8))
dev.off()
###############

### follow-up dge analysis
###############
clust.mapping <- read.csv('combined_analysis/cluster-mappings_20221104.csv')
colnames(clust.mapping) <- c('Cluster', 'Cell.Type')
sample.names <- c('KA003', 'KA005', 'KA006', 'KA007')
seurat.obj <- readRDS('combined_analysis/combined-gexp_seurat-obj.rds')
filt.counts <- readRDS('combined_analysis/combined-gexp_filt-counts.rds')
sample.name.vec <- sapply(colnames(seurat.obj), function(x) {
  strsplit(x, '\\.')[[1]][1]
})
clust.vec <- seurat.obj$seurat_clusters
## cpm normalize counts
cpm.counts <- cpm_norm(filt.counts, l2 = TRUE)
cpm.counts <- cpm.counts[, names(sample.name.vec)]
## 9+18 versus other grouped clusters
p.val.thresh <- 0.05
num.genes <- 50
test.clusts <- c(9, 18)
test.cells <- names(clust.vec)[which(clust.vec %in% test.clusts)]
ref.groups <- c('Myeloid', 'Malignant', 'Lymphoid', 'Fibroblast', 'Endothelial')
clust.9.18.dge <- list()
for (rg in ref.groups) {
  print(rg)
  # identify reference cells
  ref.clusts <- clust.mapping$Cluster[which(clust.mapping$Cell.Type == rg)]
  ref.cells <- names(clust.vec)[which(clust.vec %in% ref.clusts)]
  # run dge
  dge.results <- apply(cpm.counts, 1, function(x) {
    test.obj <- wilcox.test(x[test.cells], x[ref.cells], alternative = 'two.sided')
    rbs.cor <- 2 * test.obj$statistic / (length(test.cells) * length(ref.cells)) - 1
    return(c(test.obj$p.val, rbs.cor))
  })
  dge.results <- t(dge.results)
  dge.results <- as.data.frame(dge.results)
  colnames(dge.results) <- c('pval', 'rbsc')
  # remove na rows
  no.na.rows <- which(apply(dge.results, 1, function(x) {length(which(is.na(x))) == 0}))
  dge.results <- dge.results[no.na.rows,]
  # correct p-values
  dge.results$pval <- p.adjust(dge.results$pval, method = 'BH')  
  # add to list
  clust.9.18.dge[[rg]] <- dge.results
  # subset to sig genes
  sig.genes <- which(dge.results[,'pval'] < p.val.thresh)
  sig.dge <- dge.results[sig.genes,]
  # make pos gene df
  pos.genes <- which(sig.dge$rbsc > 0)
  if (length(pos.genes) > 0) {
    pos.df <- sig.dge[pos.genes,]
    pos.df <- pos.df[order(pos.df$rbsc, decreasing = TRUE),]
    pos.df <- pos.df[1:min(nrow(pos.df), num.genes),]
  }
  # make neg gene df
  neg.genes <- which(sig.dge$rbsc < 0)
  if (length(neg.genes) > 0) {
    neg.df <- sig.dge[neg.genes,]
    neg.df <- neg.df[order(neg.df$rbsc, decreasing = TRUE),]
    neg.df <- neg.df[1:min(nrow(neg.df), num.genes),]
  }
  # write to table
  if (length(pos.genes) > 0) {
    if (length(neg.genes) > 0) {
      combine.df <- rbind(pos.df, neg.df)
    }
  } else if (length(neg.genes > 0)) {
    combine.df <- neg.df
  } else {
    combine.df <- NA
  }
  write.csv(combine.df, file = paste('combined_analysis/cluster_followup/clust-9-18_v-', rg, '_dge.csv', sep = ''),
            row.names = TRUE, quote = FALSE)
}
saveRDS(clust.9.18.dge, file = 'combined_analysis/cluster_followup/clust-9-18_dge.rds')
# mean MT percentage for each cluster group
data("mt_genes")
mtg <- intersect(mt_genes$mur.symb, rownames(filt.counts))
mt.per <- apply(filt.counts, 2, function(x) {sum(x[mtg]) / sum(x)})
clust.mt.per <- sapply(sort(unique(clust.mapping$Cell.Type)), function(x) {
  clust.names <- clust.mapping$Cluster[which(clust.mapping$Cell.Type == x)]
  clust.samps <- names(clust.vec)[which(clust.vec %in% clust.names)]
  clust.mt <- mean(mt.per[clust.samps])
  return(clust.mt)
})
## 25 versus selected groups
###############

#######################################################################################################################################

### sample-specific gene expression clustering + plots
###############
for (sn in sample.names) {
  print(sn)
  # make directory
  dir.create(paste('sample_analysis/', sn, '/gexp_analysis/', sep = ''))
  dir.create(paste('sample_analysis/', sn, '/gexp_analysis/cluster_plots', sep = ''))
  # read filtered count matrix
  filt.counts <- readRDS(paste('sample_analysis/', sn, '/', sn, '_filt-counts.rds', sep = ''))
  # seurat analysis
  sn.seurat <- CreateSeuratObject(counts = filt.counts, assay = 'RNA', project = sn)
  sn.seurat <- PercentageFeatureSet(sn.seurat, pattern = "^mt-", col.name = "percent.mt")
  sn.seurat <- SCTransform(sn.seurat, vars.to.regress = "percent.mt", verbose = FALSE)
  sn.seurat <- RunPCA(sn.seurat, verbose = FALSE)
  sn.seurat <- RunUMAP(sn.seurat, dims = 1:30, verbose = FALSE)
  sn.seurat <- FindNeighbors(sn.seurat, dims = 1:30, verbose = FALSE)
  sn.seurat <- FindClusters(sn.seurat, verbose = FALSE)
  saveRDS(sn.seurat, file = paste('sample_analysis/', sn, '/gexp_analysis/', 
                                  sn, '_seurat.rds', sep = ''))
  # marker analysis
  cluster.markers <- list()
  for (clust.name in levels(sn.seurat$seurat_clusters)) {
    print(clust.name)
    cn.markers <- FindMarkers(sn.seurat, ident.1 = clust.name)
    cn.markers <- cn.markers[order(cn.markers$p_val_adj, -abs(cn.markers$avg_log2FC)),]
    cluster.markers[[as.character(clust.name)]] <- cn.markers
  }
  saveRDS(cluster.markers, file = paste('sample_analysis/', sn, '/gexp_analysis/', 
                                  sn, '_seurat-markers.rds', sep = ''))
  # singleR analysis
  sn.singleR <- SingleR(test = filt.counts, ref = mouse.ref,
                        labels = mouse.ref$label.main, assay.type.test = 1)
  saveRDS(sn.singleR, file = paste('sample_analysis/', sn, '/gexp_analysis/', 
                                        sn, '_singleR.rds', sep = ''))
  singleR.vec <- sn.singleR$first.labels
  names(singleR.vec) <- colnames(filt.counts)
  # umap plots
  plot.df <- data.frame('UMAP1' = sn.seurat@reductions$umap@cell.embeddings[,1],
                        'UMAP2' = sn.seurat@reductions$umap@cell.embeddings[,2],
                        'Cluster' = sn.seurat$seurat_clusters,
                        'SingleR' = singleR.vec)
  ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = Cluster), size = 0.5) +
    ggtitle(paste(sn, ': Cluster', sep = '')) +
    guides(colour = guide_legend(override.aes = list(size = 3)))
  ggsave(paste('sample_analysis/', sn, '/gexp_analysis/cluster_plots/', 
               sn, '_seurat-clust-umap.jpg', sep = ''),
         height = 8, width = 10, units = 'in', dpi = 300)
  ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = SingleR), size = 0.5) +
    ggtitle(paste(sn, ': SingleR', sep = '')) +
    guides(colour = guide_legend(override.aes = list(size = 3)))
  ggsave(paste('sample_analysis/', sn, '/gexp_analysis/cluster_plots/', 
               sn, '_singleR-umap.jpg', sep = ''),
         height = 8, width = 10, units = 'in', dpi = 300)
  # heatmap
  cell.order <- names(sort(sn.seurat$seurat_clusters))
  num.markers <- 10
  marker.set <- lapply(names(cluster.markers), function(x) {
    cm.markers <- cluster.markers[[x]]
    c.markers <- cm.markers[order(cm.markers$avg_log2FC, decreasing = TRUE),]
    return(data.frame('gene' = rownames(c.markers)[1:num.markers],
                      'group' = rep(x, num.markers)))
  })
  marker.set <- Reduce(rbind, marker.set)
  marker.set <- marker.set[which(marker.set$gene %in% rownames(sn.seurat@assays$SCT@scale.data)),]
  marker.set$group <- as.numeric(marker.set$group)
  # set colors
  clust.colors <- group_colors(length(unique(sn.seurat$seurat_clusters)))
  names(clust.colors) <- sort(unique(sn.seurat$seurat_clusters))
  singleR.colors <- group_colors(length(unique(singleR.vec)))
  names(singleR.colors) <- sort(unique(singleR.vec))
  # create annotations
  column.annot <- HeatmapAnnotation('Cluster' = sn.seurat$seurat_clusters[cell.order],
                                    'SingleR' = singleR.vec[cell.order],
                                    col = list('Cluster' = clust.colors, 'SingleR' = singleR.colors))
  marker.annot <- rowAnnotation('Cluster' = marker.set[,2],
                                col = list('Cluster' = clust.colors), 
                                show_annotation_name = FALSE,
                                show_legend = FALSE) 
  # set gaps
  col.gaps <- sn.seurat$seurat_clusters[cell.order]
  marker.gaps <- marker.set[,2]
  # plot markers
  plot.mat <- as.matrix(sn.seurat@assays$SCT@scale.data)[marker.set[,1], cell.order]
  plot.mat <- t(apply(plot.mat, 1, scale))
  col.breaks <- quantile_breaks(plot.mat)
  col.fun <- circlize::colorRamp2(col.breaks, color_levels('gexp', length(col.breaks)))
  jpeg(paste('sample_analysis/', sn, '/gexp_analysis/cluster_plots/', 
             sn, '_seurat-clust-heatmap.jpg', sep = ''),
       height = 20, width = 15, units = 'in', res = 300)
  print(Heatmap(plot.mat, name = 'GEXP',
          col = col.fun,
          top_annotation = column.annot, column_split = col.gaps,
          left_annotation = marker.annot, row_split = marker.gaps,
          cluster_rows = FALSE, cluster_columns = FALSE,
          show_row_names = TRUE, show_column_names = FALSE,
          column_title = paste(sn, ': Seurat Cluster Markers', sep = ''), row_title = NULL,
          row_names_gp = gpar(fontsize = 8)))
  dev.off()
}
###############

### plot w/ curated marker set for each sample [added-on to original sample-specific analysis]
###############
for (sn in sample.names) {
  # load objects
  filt.counts <- readRDS(paste('sample_analysis/', sn, '/', sn, '_filt-counts.rds', sep = ''))
  sn.seurat <- readRDS(paste('sample_analysis/', sn, '/gexp_analysis/', 
                             sn, '_seurat.rds', sep = ''))
  sn.singleR <- readRDS(paste('sample_analysis/', sn, '/gexp_analysis/', 
                              sn, '_singleR.rds', sep = ''))
  singleR.vec <- sn.singleR$first.labels
  names(singleR.vec) <- colnames(filt.counts)
  sct.mat <- as.matrix(sn.seurat@assays$SCT@scale.data)
  # create cpm
  cpm.counts <- cpm_norm(filt.counts, l2 = TRUE)
  # create plot data
  cell.order <- names(sort(sn.seurat$seurat_clusters))
  marker.use.inds <- which(ct.markers$Gene %in% rownames(sct.mat))
  plot.mat <- as.matrix(sct.mat[ct.markers$Gene[marker.use.inds], cell.order])
  plot.mat <- t(apply(plot.mat, 1, scale))
  # set colors
  col.breaks <- quantile_breaks(plot.mat)
  col.fun <- circlize::colorRamp2(col.breaks, color_levels('gexp', length(col.breaks)))
  clust.colors <- group_colors(length(unique(sn.seurat$seurat_clusters)))
  names(clust.colors) <- sort(unique(sn.seurat$seurat_clusters))
  singleR.colors <- group_colors(length(unique(singleR.vec)))
  names(singleR.colors) <- sort(unique(singleR.vec))
  # column annotation
  column.annot <- HeatmapAnnotation('Cluster' = sn.seurat$seurat_clusters[cell.order],
                                    'SingleR' = singleR.vec[cell.order],
                                    col = list('Cluster' = clust.colors, 'SingleR' = singleR.colors))
  col.gaps <- sn.seurat$seurat_clusters[cell.order]
  # row annotation
  row.annot <- rowAnnotation('Cell Type' = ct.markers$CT[marker.use.inds],
                             col = list('Cell Type' = ct.colors), 
                             show_annotation_name = FALSE,
                             show_legend = TRUE)
  row.gaps <- ct.markers$CT[marker.use.inds]
  # make plot
  jpeg(paste('sample_analysis/', sn, '/gexp_analysis/cluster_plots/', 
             sn, '_seurat-clust-heatmap_ct-markers.jpg', sep = ''),
       height = 10, width = 10, units = 'in', res = 300)
  print(Heatmap(plot.mat, name = 'GEXP',
                col = col.fun,
                top_annotation = column.annot, column_split = col.gaps,
                left_annotation = row.annot, row_split = row.gaps,
                cluster_rows = FALSE, cluster_columns = FALSE,
                show_row_names = TRUE, show_column_names = FALSE,
                column_title = paste(sn, ': Cell Type Markers', sep = ''), row_title = NULL,
                row_names_gp = gpar(fontsize = 8)))
  dev.off()
}
###############

### umap plots for each sample [added-on to original sample-specific analysis]
###############
for (sn in sample.names) {
  # load objects
  filt.counts <- readRDS(paste('sample_analysis/', sn, '/', sn, '_filt-counts.rds', sep = ''))
  sn.seurat <- readRDS(paste('sample_analysis/', sn, '/gexp_analysis/', 
                             sn, '_seurat.rds', sep = ''))
  sn.singleR <- readRDS(paste('sample_analysis/', sn, '/gexp_analysis/', 
                              sn, '_singleR.rds', sep = ''))
  singleR.vec <- sn.singleR$first.labels
  names(singleR.vec) <- colnames(filt.counts)
  # create plot df
  plot.df <- data.frame('UMAP1' = sn.seurat@reductions$umap@cell.embeddings[,1],
                        'UMAP2' = sn.seurat@reductions$umap@cell.embeddings[,2],
                        'Cluster' = sn.seurat$seurat_clusters,
                        'SingleR' = singleR.vec[names(sn.seurat$seurat_clusters)])
  # cluster plot
  
}
###############

#######################################################################################################################################

### gexp clustering - INCLUDES SALINE
###############
## load count matrices and combine
sample.counts <- list()
for (sn in sample.names) {
  c.mat <- readRDS(paste('sample_analysis/', sn, '/', sn, '_filt-counts.rds', sep = ''))
  sample.counts[[sn]] <- c.mat
}
shared.genes <- Reduce(intersect, lapply(sample.counts, function(x) {rownames(x)}))
counts.mat <- Reduce(cbind, lapply(sample.counts, function(x) {x[shared.genes,]}))
sample.name.vec <- rep(sample.names, sapply(sample.counts, ncol))
sample.group.vec <- rep(sample.table$Group, sapply(sample.counts, ncol))
colnames(counts.mat) <- paste(sample.name.vec, colnames(counts.mat), sep = '.')
saveRDS(counts.mat, file = 'combined_analysis/combined-gexp_filt-counts.rds')
## create seurat object and cluster
seurat.obj <- CreateSeuratObject(counts = counts.mat, assay = 'RNA', project = 'combo.seur')
seurat.obj <- PercentageFeatureSet(seurat.obj, pattern = "^mt-", col.name = "percent.mt")
seurat.obj <- SCTransform(seurat.obj, vars.to.regress = "percent.mt", verbose = FALSE)
seurat.obj <- RunPCA(seurat.obj, verbose = FALSE)
seurat.obj <- RunUMAP(seurat.obj, dims = 1:30, verbose = FALSE)
seurat.obj <- FindNeighbors(seurat.obj, dims = 1:30, verbose = FALSE)
seurat.obj <- FindClusters(seurat.obj, verbose = FALSE)
saveRDS(seurat.obj, file = 'combined_analysis/combined-gexp_seurat-obj.rds')
## seurat clust markers
cluster.markers <- list()
for (clust.name in levels(seurat.obj$seurat_clusters)) {
  print(clust.name)
  cn.markers <- FindMarkers(seurat.obj, ident.1 = clust.name)
  cn.markers <- cn.markers[order(cn.markers$p_val_adj, -abs(cn.markers$avg_log2FC)),]
  cluster.markers[[as.character(clust.name)]] <- cn.markers
}
saveRDS(cluster.markers, file = 'combined_analysis/combined-gexp_seurat-clust-markers.rds')
# save table of top and bottom markers
top.markers <- lapply(cluster.markers, function(x) {
  c.markers <- x[order(x$avg_log2FC, decreasing = TRUE),]
  return(data.frame('gene' = rownames(c.markers)[1:50],
                    'log2FC' = c.markers$avg_log2FC[1:50]))
})
top.marker.df <- Reduce(cbind, top.markers)
colnames(top.marker.df) <- paste(rep(paste('c', names(top.markers), sep = ''), each = 2),
                                 colnames(top.marker.df), sep = '.')
write.table(top.marker.df, file = 'combined_analysis/combined-gexp_seurat-clust_top-markers.csv',
            sep = ',', col.names = TRUE, row.names = FALSE, quote = FALSE)
bot.markers <- lapply(cluster.markers, function(x) {
  c.markers <- x[order(x$avg_log2FC, decreasing = FALSE),]
  return(data.frame('gene' = rownames(c.markers)[1:50],
                    'log2FC' = c.markers$avg_log2FC[1:50]))
})
bot.marker.df <- Reduce(cbind, bot.markers)
colnames(bot.marker.df) <- paste(rep(paste('c', names(bot.markers), sep = ''), each = 2),
                                 colnames(bot.marker.df), sep = '.')
write.table(bot.marker.df, file = 'combined_analysis/combined-gexp_seurat-clust_bot-markers.csv',
            sep = ',', col.names = TRUE, row.names = FALSE, quote = FALSE)
## singleR
singleR.obj <- SingleR(test = counts.mat, ref = mouse.ref,
                       labels = mouse.ref$label.main, assay.type.test = 1)
saveRDS(singleR.obj, file = 'combined_analysis/combined-gexp_singleR.rds')
singleR.vec <- singleR.obj$first.labels
names(singleR.vec) <- colnames(counts.mat)
## make heatmap
###############



