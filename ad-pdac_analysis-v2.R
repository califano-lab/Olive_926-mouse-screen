setwd('C://Users/Lukas Vlahos/OneDrive/ac_lab/ad_pdac/')
library(PISCES)
library(Seurat)
library(SingleR)
mouse.ref <- MouseRNAseqData()
sample.table <- read.table('sample-sheet.tsv', sep = '\t', header = TRUE)
sample.names <- sample.table$Sample
library(ggplot2)
library(ComplexHeatmap)
library(reshape2)
fibro.genes <- c('Gli1', 'Gli2', 'Gli3', 'Wif1', 'Ptch1')

### load general data
###############
## ct markers
ct.markers <- read.csv('pdac_ct-cluster-markers.csv', header = TRUE)
colnames(ct.markers) <- c('Gene', 'CT')
ct.colors <- group_colors(length(unique(ct.markers$CT)), offset = 60)
names(ct.colors) <- unique(ct.markers$CT)
## fibro markers
fibro.class.markers <- read.csv('pdac_fibro-markers.csv', header = TRUE) 
colnames(fibro.class.markers) <- c('Gene', 'Group')
###############

### initial QC for all samples
###############
dir.create('qc_expression')
## KA001
s.name <- 'KA001'
# load 10x outs
raw.counts <- Read10X(paste('cellranger_outs/', s.name, '_cellranger_count_outs/filtered_feature_bc_matrix/', sep = ''))
# generate qc plot
qc_plot(as.matrix(raw.counts), species = 'mur', genes = 'symb')
# filter count matrix
filt.counts <- qc_filt(as.matrix(raw.counts), min.depth = 2500, max.depth = 25000, 
                       max.genes = 5000, max.mt = 0.2, species = 'mur')
# create directory and save objects
saveRDS(filt.counts, file = paste('qc_expression/', s.name, '_filt-counts.rds', sep = ''))
jpeg(filename = paste('qc_expression/', s.name, '_qc.jpg', sep = ''),
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
                       max.genes = 6000, max.mt = 0.25, species = 'mur')
# create directory and save objects
saveRDS(filt.counts, file = paste('qc_expression/', s.name, '_filt-counts.rds', sep = ''))
jpeg(filename = paste('qc_expression/', s.name, '_qc.jpg', sep = ''),
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
                       max.genes = 6000, max.mt = 0.25, species = 'mur')
# create directory and save objects
saveRDS(filt.counts, file = paste('qc_expression/', s.name, '_filt-counts.rds', sep = ''))
jpeg(filename = paste('qc_expression/', s.name, '_qc.jpg', sep = ''),
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
                       max.genes = 6000, max.mt = 0.25, species = 'mur')
# create directory and save objects
saveRDS(filt.counts, file = paste('qc_expression/', s.name, '_filt-counts.rds', sep = ''))
jpeg(filename = paste('qc_expression/', s.name, '_qc.jpg', sep = ''),
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
                       max.genes = 6000, max.mt = 0.25, species = 'mur')
# create directory and save objects
saveRDS(filt.counts, file = paste('qc_expression/', s.name, '_filt-counts.rds', sep = ''))
jpeg(filename = paste('qc_expression/', s.name, '_qc.jpg', sep = ''),
     height = 6, width = 8, units = 'in', res = 150)
qc_plot(as.matrix(raw.counts), species = 'mur', genes = 'symb')
dev.off()
## KA008
s.name <- 'KA008'
# load 10x outs
raw.counts <- Read10X(paste('cellranger_outs/', s.name, '_cellranger_count_outs/filtered_feature_bc_matrix/', sep = ''))
# generate qc plot
qc_plot(as.matrix(raw.counts), species = 'mur', genes = 'symb')
# filter count matrix
filt.counts <- qc_filt(as.matrix(raw.counts), min.depth = 1500, max.depth = 25000, 
                       max.genes = 5000, max.mt = 0.25, species = 'mur')
# create directory and save objects
saveRDS(filt.counts, file = paste('qc_expression/', s.name, '_filt-counts.rds', sep = ''))
jpeg(filename = paste('qc_expression/', s.name, '_qc.jpg', sep = ''),
     height = 6, width = 8, units = 'in', res = 150)
qc_plot(as.matrix(raw.counts), species = 'mur', genes = 'symb')
dev.off()
## KA009
s.name <- 'KA009'
# load 10x outs
raw.counts <- Read10X(paste('cellranger_outs/', s.name, '_cellranger_count_outs/filtered_feature_bc_matrix/', sep = ''))
# generate qc plot
qc_plot(as.matrix(raw.counts), species = 'mur', genes = 'symb')
# filter count matrix
filt.counts <- qc_filt(as.matrix(raw.counts), min.depth = 1500, max.depth = 35000, 
                       max.genes = 7500, max.mt = 0.25, species = 'mur')
# create directory and save objects
saveRDS(filt.counts, file = paste('qc_expression/', s.name, '_filt-counts.rds', sep = ''))
jpeg(filename = paste('qc_expression/', s.name, '_qc.jpg', sep = ''),
     height = 6, width = 8, units = 'in', res = 150)
qc_plot(as.matrix(raw.counts), species = 'mur', genes = 'symb')
dev.off()
###############

#############################################
## COMBINED ANALYSIS w/ NO ANCHORING ##

### combined clustering analysis
###############
sample.names <- c('KA003', 'KA005', 'KA006', 'KA007', 'KA008', 'KA009')
sample.inds <- match(sample.names, sample.table$Sample)
## load count matrices and combine
sample.counts <- list()
for (sn in sample.names) {
  print(sn)
  c.mat <- readRDS(paste('qc_expression/', sn, '_filt-counts.rds', sep = ''))
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
## pca plot to check for batch effects
sample.name.vec <- sapply(colnames(counts.mat), function(x) {strsplit(x, '\\.')[[1]][1]} )
group.name.vec <- sample.table$Group[match(sample.name.vec, sample.table$Sample)]
names(group.name.vec) <- colnames(counts.mat)
plot.df <- data.frame('PC1' = seurat.obj@reductions$pca@cell.embeddings[,1],
                      'PC2' = seurat.obj@reductions$pca@cell.embeddings[,2],
                      'Sample' = sample.name.vec,
                      'Group' = group.name.vec)
ggplot(plot.df, aes(PC1, PC2)) + geom_point(aes(color = Sample), size = 0.5) +
  labs(title = 'Combined Counts - Seurat PCA') +
  guides(colour = guide_legend(override.aes = list(size=3)))
ggsave('combined_analysis/combined-gexp_seurat-pca.jpg', height = 8, width = 10, units = 'in', dpi = 300)
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
write.table(top.marker.df, file = 'combined_analysis/combined-gexp_seurat-clust_top-markers_v2.csv',
            sep = ',', col.names = TRUE, row.names = FALSE, quote = FALSE)
bot.markers <- lapply(cluster.markers, function(x) {
  c.markers <- x[order(x$avg_log2FC, decreasing = FALSE),]
  return(data.frame('gene' = rownames(c.markers)[1:50],
                    'log2FC' = c.markers$avg_log2FC[1:50]))
})
bot.marker.df <- Reduce(cbind, bot.markers)
colnames(bot.marker.df) <- paste(rep(paste('c', names(bot.markers), sep = ''), each = 2),
                                 colnames(bot.marker.df), sep = '.')
write.table(bot.marker.df, file = 'combined_analysis/combined-gexp_seurat-clust_bot-markers_v2.csv',
            sep = ',', col.names = TRUE, row.names = FALSE, quote = FALSE)
## singleR
singleR.obj <- SingleR(test = counts.mat, ref = mouse.ref,
                       labels = mouse.ref$label.main, assay.type.test = 1)
saveRDS(singleR.obj, file = 'combined_analysis/combined-gexp_singleR.rds')
singleR.vec <- singleR.obj$first.labels
names(singleR.vec) <- colnames(counts.mat)
## get sample and group vectors
sample.name.vec <- sapply(colnames(counts.mat), function(x) {strsplit(x, '\\.')[[1]][1]} )
group.name.vec <- sample.table$Group[match(sample.name.vec, sample.table$Sample)]
names(group.name.vec) <- colnames(counts.mat)
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
ggsave('combined_analysis/cluster_plots/combined-gexp_sample-umap_v2.jpg',
       height = 8, width = 10, units = 'in', dpi = 300)
ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = Group), size = 0.5) +
  ggtitle('Combined Analysis - Group') +
  guides(colour = guide_legend(override.aes = list(size = 3)))
ggsave('combined_analysis/cluster_plots/combined-gexp_group-umap_v2.jpg',
       height = 8, width = 10, units = 'in', dpi = 300)
ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = Cluster), size = 0.5) +
  ggtitle('Combined Analysis - Cluster') +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  annotate("text", x = cluster.umap.centers$UMAP_1, y = cluster.umap.centers$UMAP_2,
           label = cluster.umap.centers$Cluster)
ggsave('combined_analysis/cluster_plots/combined-gexp_cluster-umap_v2.jpg',
       height = 8, width = 10, units = 'in', dpi = 300)
ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = SingleR), size = 0.5) +
  ggtitle('Combined Analysis - SingleR') +
  guides(colour = guide_legend(override.aes = list(size = 3)))
ggsave('combined_analysis/cluster_plots/combined-gexp_singler-umap_v2.jpg',
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
jpeg('combined_analysis/cluster_plots/combined-gexp_cluster-markers_v2.jpg',
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

### sample specific inferCNV based on combined analysis clustering
## fibroblast reference
###############
library(infercnv)
sample.names <- c('KA003', 'KA005', 'KA006', 'KA007')
filt.counts <- readRDS('combined_analysis/combined-gexp_filt-counts.rds')
seurat.obj <- readRDS('combined_analysis/combined-gexp_seurat-obj.rds')
cluster.mappings <- read.csv('combined_analysis/cluster-mapping_20221107.csv')
colnames(cluster.mappings) <- c('Cluster', 'Cell.Type', 'Sub.Type')
sample.name.vec <- sapply(colnames(seurat.obj), function(x) {
  strsplit(x, '\\.')[[1]][1]
})
## create labeling vector
clust.vec <- seurat.obj$seurat_clusters
ref.clusters <- cluster.mappings$Cluster[which(cluster.mappings$Cell.Type == 'Fibroblast')]
ref.clusters <- as.character(ref.clusters)
test.clusters <- cluster.mappings$Cluster[which(cluster.mappings$Cell.Type %in% c('Epithelial'))]
test.clusters <- as.character(test.clusters)
## perform analysis for each sample
for (sn in sample.names) {
  print(sn)
  # create directory
  out.dir <- paste('combined_analysis/infercnv/', sn, '_fibro-ref/', sep = '')
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
library(infercnv)
sample.names <- c('KA003', 'KA005', 'KA006', 'KA007')
filt.counts <- readRDS('combined_analysis/combined-gexp_filt-counts.rds')
seurat.obj <- readRDS('combined_analysis/combined-gexp_seurat-obj.rds')
cluster.mappings <- read.csv('combined_analysis/cluster-mapping_20221107.csv')
colnames(cluster.mappings) <- c('Cluster', 'Cell.Type', 'Sub.Type')
sample.name.vec <- sapply(colnames(seurat.obj), function(x) {
  strsplit(x, '\\.')[[1]][1]
})
## create labeling vector
clust.vec <- seurat.obj$seurat_clusters
ref.clusters <- cluster.mappings$Cluster[which(cluster.mappings$Cell.Type %in% c('Fibroblast', 'Endothelial'))]
ref.clusters <- as.character(ref.clusters)
test.clusters <- cluster.mappings$Cluster[which(cluster.mappings$Cell.Type %in% c('Epithelial'))]
test.clusters <- as.character(test.clusters)
## perform analysis for each sample
for (sn in sample.names) {
  print(sn)
  # create directory
  out.dir <- paste('combined_analysis/infercnv/', sn, '_fibro-endo-ref/', sep = '')
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

### group frequency analysis
###############
sample.names <- c('KA003', 'KA005', 'KA006', 'KA007')
filt.counts <- readRDS('combined_analysis/combined-gexp_filt-counts.rds')
seurat.obj <- readRDS('combined_analysis/combined-gexp_seurat-obj.rds')
cluster.mappings <- read.csv('combined_analysis/cluster-mapping_20221107.csv')
colnames(cluster.mappings) <- c('Cluster', 'Cell.Type', 'Sub.Type')
## create label vectors
sample.name.vec <- sapply(colnames(seurat.obj), function(x) {strsplit(x, '\\.')[[1]][1]} )
group.name.vec <- sample.table$Group[match(sample.name.vec, sample.table$Sample)]
names(group.name.vec) <- names(sample.name.vec)
clust.vec <- seurat.obj$seurat_clusters
ct.vec <- plyr::mapvalues(clust.vec, from = cluster.mappings$Cluster, to = cluster.mappings$Cell.Type)
## write ct confusion matrices for later reference
sample.ct.table <- table(sample.name.vec, ct.vec[names(sample.name.vec)])
write.csv(sample.ct.table, 'combined_analysis/ct-frequency_analysis/sample-ct_table.csv')
group.ct.table <- table(group.name.vec, ct.vec[names(group.name.vec)])
write.csv(group.ct.table, 'combined_analysis/ct-frequency_analysis/group-ct_table.csv')
## comparison of vehicle treatments
test.samps <- names(which(group.name.vec == 'Vehicle'))
test.table <- table(sample.name.vec[test.samps], ct.vec[test.samps])
test.obj <- chisq.test(x = test.table)
test.adjp <- matrix(p.adjust(method = 'BH', p = pnorm(abs(test.obj$stdres), lower.tail = FALSE)), nrow = 2)
test.adjp <- test.adjp[1,]
test.sign <- sign(test.obj$stdres)
test.list <- list('adj.p' = test.adjp, 'sign' = test.sign)
test.df <- rbind(test.sign, test.adjp)
rownames(test.df) <- c(paste(rownames(test.sign), 'sign', sep = '.'), 'P-val')
saveRDS(test.list, file = 'combined_analysis/ct-frequency_analysis/vehicle-similarity_chisq.rds')
write.csv(test.df, file = 'combined_analysis/ct-frequency_analysis/vehicle-similarity_chisq.csv')
## comaprison of drug treatments
test.samps <- names(which(group.name.vec == '926'))
test.table <- table(sample.name.vec[test.samps], ct.vec[test.samps])
test.obj <- chisq.test(x = test.table)
test.adjp <- matrix(p.adjust(method = 'BH', p = pnorm(abs(test.obj$stdres), lower.tail = FALSE)), nrow = 2)
test.adjp <- test.adjp[1,]
test.sign <- sign(test.obj$stdres)
test.list <- list('adj.p' = test.adjp, 'sign' = test.sign)
test.df <- rbind(test.sign, test.adjp)
rownames(test.df) <- c(paste(rownames(test.sign), 'sign', sep = '.'), 'P-val')
saveRDS(test.list, file = 'combined_analysis/ct-frequency_analysis/926-similarity_chisq.rds')
write.csv(test.df, file = 'combined_analysis/ct-frequency_analysis/926-similarity_chisq.csv')
## comparison between vehicle and drug
test.samps <- names(group.name.vec)
test.table <- table(group.name.vec[test.samps], ct.vec[test.samps])
test.obj <- chisq.test(x = test.table)
test.adjp <- matrix(p.adjust(method = 'BH', p = pnorm(abs(test.obj$stdres), lower.tail = FALSE)), nrow = 2)
test.adjp <- test.adjp[1,]
test.sign <- sign(test.obj$stdres)
test.list <- list('adj.p' = test.adjp, 'sign' = test.sign)
test.df <- rbind(test.sign, test.adjp)
rownames(test.df) <- c(paste(rownames(test.sign), 'sign', sep = '.'), 'P-val')
saveRDS(test.list, file = 'combined_analysis/ct-frequency_analysis/926-v-vehicle_chisq.rds')
write.csv(test.df, file = 'combined_analysis/ct-frequency_analysis/926-v-vehicle_chisq.csv')
## pairwise test approach
treated.samps <- intersect(sample.names, sample.table$Sample[which(sample.table$Group == '926')])
vehicle.samps <- intersect(sample.names, sample.table$Sample[which(sample.table$Group == 'Vehicle')])
test.obj.list <- list()
for (ts in treated.samps) {
  for (vs in vehicle.samps) {
    # get samples
    test.samps <- names(sample.name.vec)[which(sample.name.vec %in% c(ts, vs))]
    # run test
    test.table <- table(sample.name.vec[test.samps], ct.vec[test.samps])
    print(test.table)
    test.obj <- chisq.test(x = test.table)
    test.adjp <- matrix(p.adjust(method = 'BH', p = pnorm(abs(test.obj$stdres), lower.tail = FALSE)), nrow = 2)
    test.adjp <- test.adjp[1,]
    test.sign <- sign(test.obj$stdres)
    test.list <- list('adj.p' = test.adjp, 'sign' = test.sign)
    test.df <- rbind(test.sign, test.adjp)
    rownames(test.df) <- c(paste(rownames(test.sign), 'sign', sep = '.'), 'P-val')
    # add to list
    test.obj.list[[paste(ts, 'v', vs, sep = '.')]] <- list('test.list' = test.list, 'test.df' = test.df)
  }
}
# integrate into final df
integrated.sign.df <- Reduce('+', lapply(test.obj.list, function(x) {x$test.df[1:2,]})) / 4
integrated.pval.df <- Reduce('+', lapply(test.obj.list, function(x) {x$test.df[3,]})) / sqrt(4)
integrated.df <- rbind(integrated.sign.df, integrated.pval.df)
rownames(integrated.df) <- c('926.sign', 'Vehicle.sign', 'PVal')
write.csv(integrated.df, file = 'combined_analysis/ct-frequency_analysis/pairwise_926-v-vehicle_chisq.csv')
# save as RDS
test.save.obj <- list('test.objects' = test.obj.list, 'integrated.df' = integrated.df)
saveRDS(test.obj.list, file = 'combined_analysis/ct-frequency_analysis/pairwise_926-v-vehicle_chisq.rds')
###############

#############################################





#############################################
## COMBINED ANALYSIS w/ ANCHORING ##
## INCLUDES KA008 and KA009 ##

### combined analysis w/ anchoring
###############
sample.names <- c('KA003', 'KA005', 'KA006', 'KA007', 'KA008', 'KA009')
sample.inds <- match(sample.names, sample.table$Sample)
## make directory
dir.create('anchor-combined_analysis/')
## load count matrices and combine
seurat.obj.list <- list()
for (sn in sample.names) {
  print(sn)
  # load counts
  c.mat <- readRDS(paste('qc_expression/', sn, '_filt-counts.rds', sep = ''))
  colnames(c.mat) <- paste(sn, colnames(c.mat), sep = '.')
  # create seurat object
  sn.seurat.obj <- CreateSeuratObject(counts = c.mat, assay = 'RNA', project = sn)
  sn.seurat.obj <- PercentageFeatureSet(sn.seurat.obj, pattern = "^mt-", col.name = "percent.mt")
  sn.seurat.obj <- SCTransform(sn.seurat.obj, vars.to.regress = "percent.mt", verbose = FALSE)
  print(head(colnames(sn.seurat.obj)))
  # add to list
  seurat.obj.list[[sn]] <- sn.seurat.obj
}
saveRDS(seurat.obj.list, file = 'anchor-combined_analysis/seurat-obj_list.rds')
## anchor
integration.features<- SelectIntegrationFeatures(object.list = seurat.obj.list, nfeatures = 4000)
integration.anchors <- FindIntegrationAnchors(object.list = seurat.obj.list, anchor.features = integration.features)
int.seurat.obj <- IntegrateData(anchorset = integration.anchors)
saveRDS(int.seurat.obj, file = 'anchor-combined_analysis/anchored_seurat-obj.rds')
remove(seurat.obj.list)
## clustering analysis
int.seurat.obj <- SCTransform(int.seurat.obj, vars.to.regress = "percent.mt", verbose = FALSE)
int.seurat.obj <- RunPCA(int.seurat.obj, verbose = FALSE)
int.seurat.obj <- RunUMAP(int.seurat.obj, dims = 1:30, verbose = FALSE)
int.seurat.obj <- FindNeighbors(int.seurat.obj, dims = 1:30, verbose = FALSE)
int.seurat.obj <- FindClusters(int.seurat.obj, verbose = FALSE)
saveRDS(int.seurat.obj, file = 'anchor-combined_analysis/anchored_seurat-obj.rds')
## cluster markers
cluster.markers <- list()
for (clust.name in levels(int.seurat.obj$seurat_clusters)) {
  print(clust.name)
  cn.markers <- FindMarkers(int.seurat.obj, ident.1 = clust.name)
  cn.markers <- cn.markers[order(cn.markers$p_val_adj, -abs(cn.markers$avg_log2FC)),]
  cluster.markers[[as.character(clust.name)]] <- cn.markers
}
saveRDS(cluster.markers, file = 'anchor-combined_analysis/anchored_seurat-clust-markers.rds')
# save table of top and bottom markers
top.markers <- lapply(cluster.markers, function(x) {
  c.markers <- x[order(x$avg_log2FC, decreasing = TRUE),]
  return(data.frame('gene' = rownames(c.markers)[1:50],
                    'log2FC' = c.markers$avg_log2FC[1:50]))
})
top.marker.df <- Reduce(cbind, top.markers)
colnames(top.marker.df) <- paste(rep(paste('c', names(top.markers), sep = ''), each = 2),
                                 colnames(top.marker.df), sep = '.')
write.table(top.marker.df, file = 'anchor-combined_analysis/combined-gexp_seurat-clust_top-markers_v2.csv',
            sep = ',', col.names = TRUE, row.names = FALSE, quote = FALSE)
bot.markers <- lapply(cluster.markers, function(x) {
  c.markers <- x[order(x$avg_log2FC, decreasing = FALSE),]
  return(data.frame('gene' = rownames(c.markers)[1:50],
                    'log2FC' = c.markers$avg_log2FC[1:50]))
})
bot.marker.df <- Reduce(cbind, bot.markers)
colnames(bot.marker.df) <- paste(rep(paste('c', names(bot.markers), sep = ''), each = 2),
                                 colnames(bot.marker.df), sep = '.')
write.table(bot.marker.df, file = 'anchor-combined_analysis/combined-gexp_seurat-clust_bot-markers_v2.csv',
            sep = ',', col.names = TRUE, row.names = FALSE, quote = FALSE)
## singleR
counts.mat <- readRDS(file = 'combined_analysis/combined-gexp_filt-counts.rds')
singleR.obj <- readRDS('combined_analysis/combined-gexp_singleR.rds')
singleR.vec <- singleR.obj$first.labels
names(singleR.vec) <- colnames(counts.mat)
## get sample and group vectors
sample.name.vec <- sapply(colnames(counts.mat), function(x) {strsplit(x, '\\.')[[1]][1]} )
group.name.vec <- sample.table$Group[match(sample.name.vec, sample.table$Sample)]
names(group.name.vec) <- colnames(counts.mat)
rm(counts.mat)
## umap plots
plot.df <- data.frame('UMAP1' = int.seurat.obj@reductions$umap@cell.embeddings[,1],
                      'UMAP2' = int.seurat.obj@reductions$umap@cell.embeddings[,2],
                      'Sample' = sample.name.vec, 'Group' = group.name.vec,
                      'Cluster' = int.seurat.obj$seurat_clusters,
                      'SingleR' = singleR.vec)
cluster.umap.centers <- lapply(sort(unique(int.seurat.obj$seurat_clusters)), function(x) {
  clust.samps <- which(int.seurat.obj$seurat_clusters == x)
  mean.vals <- colMeans(int.seurat.obj@reductions$umap@cell.embeddings[clust.samps, 1:2])
  return(mean.vals)
})
cluster.umap.centers <- Reduce(rbind, cluster.umap.centers)
cluster.umap.centers <- as.data.frame(cluster.umap.centers)
cluster.umap.centers$Cluster <- sort(unique(int.seurat.obj$seurat_clusters))
# sample plot
ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = Sample), size = 0.5) +
  ggtitle('Anchored Combined Analysis - Sample') +
  guides(colour = guide_legend(override.aes = list(size = 3)))
ggsave('anchor-combined_analysis/cluster_plots/combined-gexp_sample-umap_v2.jpg',
       height = 8, width = 10, units = 'in', dpi = 300)
ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = Group), size = 0.5) +
  ggtitle('Anchored Combined Analysis - Group') +
  guides(colour = guide_legend(override.aes = list(size = 3)))
ggsave('anchor-combined_analysis/cluster_plots/combined-gexp_group-umap_v2.jpg',
       height = 8, width = 10, units = 'in', dpi = 300)
ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = Cluster), size = 0.5) +
  ggtitle('Anchored Combined Analysis - Cluster') +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  annotate("text", x = cluster.umap.centers$UMAP_1, y = cluster.umap.centers$UMAP_2,
           label = cluster.umap.centers$Cluster)
ggsave('anchor-combined_analysis/cluster_plots/combined-gexp_cluster-umap_v2.jpg',
       height = 8, width = 10, units = 'in', dpi = 300)
ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = SingleR), size = 0.5) +
  ggtitle('Anchored Combined Analysis - SingleR') +
  guides(colour = guide_legend(override.aes = list(size = 3)))
ggsave('anchor-combined_analysis/cluster_plots/combined-gexp_singler-umap_v2.jpg',
       height = 8, width = 10, units = 'in', dpi = 300)
## map names
sct.mat <- as.matrix(int.seurat.obj@assays$SCT@scale.data)
colnames(sct.mat) <- colnames(counts.mat)
clust.vec <- int.seurat.obj$seurat_clusters
names(clust.vec) <- colnames(counts.mat)
## marker heatmap
cell.order <- names(sort(clust.vec))
marker.use.inds <- which(ct.markers$Gene %in% rownames(sct.mat))
plot.mat <- sct.mat[ct.markers$Gene[marker.use.inds], cell.order]
plot.mat <- t(apply(plot.mat, 1, scale))
# set colors
col.breaks <- quantile_breaks(plot.mat)
col.fun <- circlize::colorRamp2(col.breaks, color_levels('gexp', length(col.breaks)))
clust.colors <- group_colors(length(unique(int.seurat.obj$seurat_clusters)))
names(clust.colors) <- sort(unique(int.seurat.obj$seurat_clusters))
singleR.colors <- group_colors(length(unique(singleR.vec)))
names(singleR.colors) <- sort(unique(singleR.vec))
sample.colors <- group_colors(length(unique(sample.name.vec)))
names(sample.colors) <- sort(unique(sample.name.vec))
group.colors <- c('darkgrey', 'lightgrey')
names(group.colors) <- sort(unique(group.name.vec))
# column annotation
column.annot <- HeatmapAnnotation('Cluster' = clust.vec[cell.order],
                                  'SingleR' = singleR.vec[cell.order],
                                  'Sample' = sample.name.vec[cell.order],
                                  'Group' = group.name.vec[cell.order],
                                  col = list('Cluster' = clust.colors, 'SingleR' = singleR.colors,
                                             'Sample' = sample.colors, 'Group' = group.colors))
col.gaps <- clust.vec[cell.order]
# row annotation
row.annot <- rowAnnotation('Cell Type' = ct.markers$CT[marker.use.inds],
                           col = list('Cell Type' = ct.colors), 
                           show_annotation_name = FALSE,
                           show_legend = TRUE)
row.gaps <- ct.markers$CT[marker.use.inds]
# make plot
jpeg('anchor-combined_analysis/cluster_plots/combined-gexp_cluster-markers_v2.jpg',
     height = 10, width = 15, units = 'in', res = 300)
print(Heatmap(plot.mat, name = 'GEXP',
              col = col.fun,
              top_annotation = column.annot, column_split = col.gaps,
              left_annotation = row.annot, row_split = row.gaps,
              cluster_rows = FALSE, cluster_columns = FALSE,
              show_row_names = TRUE, show_column_names = FALSE,
              column_title = 'Anchored Combined Analysis - Supervised Markers', row_title = NULL,
              row_names_gp = gpar(fontsize = 8)))
dev.off()
###############

### sample specific inferCNV based on combined analysis clustering w/ anchoring
## fibroblast reference
###############
library(infercnv)
sample.names <- c('KA003', 'KA005', 'KA006', 'KA007', 'KA008', 'KA009')
int.seurat.obj <- readRDS('anchor-combined_analysis/anchored_seurat-obj.rds')
cluster.mappings <- read.csv('combined_analysis/cluster-mapping_20230124.csv')
colnames(cluster.mappings) <- c('Cluster', 'Cell.Type', 'Sub.Type')
## create labeling vector
clust.vec <- int.seurat.obj$seurat_clusters
ref.clusters <- cluster.mappings$Cluster[which(cluster.mappings$Cell.Type == 'Fibroblast')]
ref.clusters <- as.character(ref.clusters)
ref.cluster.cells <- names(clust.vec)[which(clust.vec %in% ref.clusters)]
test.clusters <- cluster.mappings$Cluster[which(cluster.mappings$Cell.Type %in% c('Epithelial'))]
test.clusters <- as.character(test.clusters)
test.cluster.cells <- names(clust.vec)[which(clust.vec %in% test.clusters)]
## perform analysis for each sample
for (sn in sample.names) {
  print(sn)
  # create directory
  out.dir <- paste('anchor-combined_analysis/infercnv/', sn, '_fibro-ref/', sep = '')
  dir.create(out.dir)
  # load counts
  sn.counts <- readRDS(paste('qc_expression/', sn, '_filt-counts.rds', sep = ''))
  colnames(sn.counts) <- paste(sn, colnames(sn.counts), sep = '.')
  # subset to the reference and test clusters with enough data for this sn
  ref.cells <- intersect(colnames(sn.counts), ref.cluster.cells)
  test.cells <- intersect(colnames(sn.counts), test.cluster.cells)
  sn.ref.clusts <- table(clust.vec[ref.cells])
  sn.ref.clusts <- names(sn.ref.clusts)[which(sn.ref.clusts > 1)]
  sn.test.clusts <- table(clust.vec[test.cells])
  sn.test.clusts <- names(sn.test.clusts)[which(sn.test.clusts > 1)]
  # create icnv objects
  ref.cells <- intersect(colnames(sn.counts), names(clust.vec)[which(clust.vec %in% sn.ref.clusts)])
  test.cells <- intersect(colnames(sn.counts), names(clust.vec)[which(clust.vec %in% sn.test.clusts)])
  icnv.mat <- sn.counts[, c(ref.cells, test.cells)]
  icnv.ref <- clust.vec[c(ref.cells, test.cells)]
  # save objects
  write.table(icnv.mat, sep = '\t', quote = FALSE,
              file = paste(out.dir, 'combined_icnv-in.txt', sep = ''))
  write.table(icnv.ref, sep = '\t', quote = FALSE, col.names = FALSE,
              file = paste(out.dir, 'combined_icnv-annot.txt', sep = ''))
  # run infercnv
  icnv.obj <- CreateInfercnvObject(raw_counts_matrix = paste(out.dir, 'combined_icnv-in.txt', sep = ''),
                                   annotations_file = paste(out.dir, 'combined_icnv-annot.txt', sep = ''),
                                   gene_order_file = '../reference_genomes/mouse_chr_gene_locations.txt',
                                   delim = '\t', ref_group_names = sn.ref.clusts)
  icnv.cells <- run(icnv.obj, cutoff = 0.1, out_dir = out.dir,
                    cluster_by_groups = TRUE, k_obs_groups = 1, analysis_mode = 'samples', 
                    denoise = TRUE, noise_logistic = TRUE,
                    scale_data = FALSE,
                    HMM = FALSE, HMM_type = 'i6', 
                    no_plot = FALSE, plot_steps = TRUE, debug = FALSE, num_threads = 11, up_to_step = 50)
}
###############

### make combined count matrix
###############
sample.names <- c('KA003', 'KA005', 'KA006', 'KA007', 'KA008', 'KA009')
sample.inds <- match(sample.names, sample.table$Sample)
## load count matrices
sample.counts <- list()
for (sn in sample.names) {
  print(sn)
  c.mat <- readRDS(paste('qc_expression/', sn, '_filt-counts.rds', sep = ''))
  sample.counts[[sn]] <- c.mat
}
## get union genes
union.genes <- unique(unlist(sapply(sample.counts, rownames)))
## pad each matrix with missing genes
pad.sample.counts <- lapply(sample.counts, function(x) {
  missing.genes <- setdiff(union.genes, rownames(x))
  pad.mat <- matrix(0L, nrow = length(missing.genes), ncol = ncol(x))
  rownames(pad.mat) <- missing.genes
  x.mat <- rbind(x, pad.mat)
  x.mat <- x.mat[union.genes,]
  return(x.mat)
})
## modify column names
sample.name.pref <- rep(sample.names, sapply(pad.sample.counts, ncol))
## combine into final matrix
pad.counts.mat <- Reduce(cbind, pad.sample.counts)
colnames(pad.counts.mat) <- paste(sample.name.pref, colnames(pad.counts.mat), sep = '.')
saveRDS(pad.counts.mat, file = 'anchor-combined_analysis/combined-gexp_pad-filt-counts.rds')
###############

### group vehicle networks
###############
sample.names <- c('KA003', 'KA005', 'KA006', 'KA007', 'KA008', 'KA009')
filt.counts <- readRDS('anchor-combined_analysis/combined-gexp_pad-filt-counts.rds')
int.seurat.obj <- readRDS('anchor-combined_analysis/anchored_seurat-obj.rds')
cluster.mappings <- read.csv('combined_analysis/cluster-mapping_20230124.csv')
colnames(cluster.mappings) <- c('Cluster', 'Cell.Type', 'Sub.Type')
## create label vectors
sample.name.vec <- sapply(colnames(int.seurat.obj), function(x) {strsplit(x, '\\.')[[1]][1]} )
group.name.vec <- sample.table$Group[match(sample.name.vec, sample.table$Sample)]
names(group.name.vec) <- names(sample.name.vec)
clust.vec <- int.seurat.obj$seurat_clusters
ct.vec <- plyr::mapvalues(clust.vec, from = cluster.mappings$Cluster, to = cluster.mappings$Cell.Type)
## table of vehicle samps and cell types
vehicle.samps <- names(group.name.vec)[which(group.name.vec == 'Vehicle')]
vehicle.ct.table <- table(sample.name.vec[vehicle.samps], ct.vec[vehicle.samps])
## make aracne matrices
dir.create('anchor-combined_analysis/a3_mats/')
set.seed(343)
for (ct in unique(ct.vec)) {
  print(ct)
  # make matrix
  ct.samps <- intersect(vehicle.samps, names(ct.vec)[which(ct.vec == ct)])
  a3.mat <- filt.counts[, ct.samps]
  if (ncol(a3.mat) > 1000) {a3.mat <- a3.mat[, sample(ct.samps, 1000)]}
  a3.mat <- a3.mat[which(rowSums(a3.mat) > 0),]
  a3.mat <- cpm_norm(a3.mat, l2 = TRUE)
  print(dim(a3.mat))
  print(mean(colSums(a3.mat)))
  # save
  ct.name <- tolower(ct)
  saveRDS(a3.mat, paste('anchor-combined_analysis/a3_mats/', ct.name, '_a3-mat.rds', sep = ''))
}
###############

### create network lists
###############
network.names <- c('endo', 'epi', 'fibro', 'lymph', 'myeloid')
dir.create('anchor-combined_analysis/a3_nets/net-lists/')
network.mdata <- list()
for (nn in network.names) {
  print(nn)
  nn.tf <- readRDS(paste('anchor-combined_analysis/a3_nets/ad-', nn, '_tf_regulon-list.rds', sep = ''))
  nn.sig <- readRDS(paste('anchor-combined_analysis/a3_nets/ad-', nn, '_sig_regulon-list.rds', sep = ''))
  nn.surf <- readRDS(paste('anchor-combined_analysis/a3_nets/ad-', nn, '_surf_regulon-list.rds', sep = ''))
  nn.net <- c(nn.tf, nn.sig, nn.surf)
  saveRDS(nn.net, file = paste('anchor-combined_analysis/a3_nets/net-lists/ad-', nn, '_net-list.rds', sep = ''))
  mdata.vec <- c(length(nn.tf), length(nn.sig), length(nn.surf), length(nn.net))
  network.mdata[[nn]] <- mdata.vec
}
mdata.df <- as.data.frame(Reduce(rbind, network.mdata))
rownames(mdata.df) <- network.names
colnames(mdata.df) <- c('TF', 'Sig', 'Surf', 'Total')
write.csv(mdata.df, file = 'anchor-combined_analysis/a3_nets/net-size_metadata.csv', quote = FALSE)
###############

### network matched narnea analysis
###############
filt.counts <- readRDS('anchor-combined_analysis/combined-gexp_pad-filt-counts.rds')
int.seurat.obj <- readRDS('anchor-combined_analysis/anchored_seurat-obj.rds')
cluster.mappings <- read.csv('combined_analysis/cluster-mapping_20230124.csv')
colnames(cluster.mappings) <- c('Cluster', 'Cell.Type', 'Sub.Type')
## create label vectors
sample.name.vec <- sapply(colnames(int.seurat.obj), function(x) {strsplit(x, '\\.')[[1]][1]} )
group.name.vec <- sample.table$Group[match(sample.name.vec, sample.table$Sample)]
names(group.name.vec) <- names(sample.name.vec)
clust.vec <- int.seurat.obj$seurat_clusters
ct.vec <- plyr::mapvalues(clust.vec, from = cluster.mappings$Cluster, to = cluster.mappings$Cell.Type)
## mapping from ct to network names
ct.network.names <- c('myeloid', 'fibro', 'epi', 'lymph', 'endo')
names(ct.network.names) <- sort(unique(ct.vec))
## narnea analysis for each cell type - 926 vs vehicle
for (ct.name in sort(unique(ct.vec))[2:5]) {
  print(ct.name)
  ct.samps <- names(ct.vec)[which(ct.vec == ct.name)]
  dir.create(paste('anchor-combined_analysis/', ct.name, sep = ''))
  
  # create anchored object for this cell type
  cat('Anchoring...\n')
  seur.obj.list <- list()
  for (s.name in sort(unique(sample.name.vec))) {
    print(s.name)
    s.samps <- names(sample.name.vec)[which(sample.name.vec == s.name)]
    s.ct.samps <- intersect(ct.samps, s.samps)
    s.ct.counts <- filt.counts[, s.ct.samps]
    s.ct.seurat.obj <- CreateSeuratObject(counts = s.ct.counts, assay = 'RNA', project = paste(s.name, ct.name, sep = '.'))
    s.ct.seurat.obj <- PercentageFeatureSet(s.ct.seurat.obj, pattern = "^mt-", col.name = "percent.mt")
    s.ct.seurat.obj <- SCTransform(s.ct.seurat.obj, vars.to.regress = "percent.mt", verbose = FALSE)
    seur.obj.list[[s.name]] <- s.ct.seurat.obj
  }
  int.features<- SelectIntegrationFeatures(object.list = seur.obj.list, nfeatures = 4000)
  int.anchors <- FindIntegrationAnchors(object.list = seur.obj.list, anchor.features = int.features)
  s.ct.int.seurat <- IntegrateData(anchorset = int.anchors, k.weight = min(100, min(sapply(seur.obj.list, ncol))))
  s.ct.int.seurat <- SCTransform(s.ct.int.seurat, vars.to.regress = "percent.mt", verbose = FALSE)
  saveRDS(s.ct.int.seurat, file = paste('anchor-combined_analysis/', ct.name, '/',
                                        ct.name, '_int-seurat-obj.rds', sep = ''))
  rm(seur.obj.list)
  rm(s.ct.seurat.obj)
  
  # create GES
  cat('Generating GES...\n')
  sct.mat <- s.ct.int.seurat@assays$SCT@scale.data
  test.cells <- intersect(colnames(sct.mat), names(group.name.vec)[which(group.name.vec == '926')])
  ref.cells <- intersect(colnames(sct.mat), names(group.name.vec)[which(group.name.vec == 'Vehicle')])
  test.mat <- sct.mat[, test.cells]
  ref.mat <- sct.mat[, ref.cells]
  vs.vehicle.ges <- cbind(scale_ges(test.mat, ref.mat), t(scale(t(ref.mat))))
  saveRDS(vs.vehicle.ges, file = paste('anchor-combined_analysis/', ct.name, '/',
                                       ct.name, '_int-seurat_vs-vehicle-ges.rds', sep = ''))
  rm(test.mat)
  rm(ref.mat)
  rm(sct.mat)
  
  # run narnea
  cat('Running NaRnEA...\n')
  ct.net.list <- readRDS(paste('anchor-combined_analysis/a3_nets/net-lists/ad-', 
                               ct.network.names[ct.name], '_net-list.rds', sep = ''))
  narnea.obj <- matrix_narnea(vs.vehicle.ges, ct.net.list)
  has.na.row <- apply(narnea.obj$PES, 1, function(x) {length(which(is.na(x)) > 0)})
  narnea.obj <- list('PES' = narnea.obj$PES[which(has.na.row == 0),],
                     'NES' = narnea.obj$NES[which(has.na.row == 0),])
  saveRDS(narnea.obj, file = paste('anchor-combined_analysis/', ct.name, '/',
                                   ct.name, '_int-seurat_vs-vehicle-ges_narnea.rds', sep = ''))
  
  # narnea analysis
  cat('NaRnEA analysis...\n')
  narnea.dist <- dist(t(narnea.obj$PES))
  narnea.mds <- make_mds(narnea.dist)
  narnea.clust <- louvain_k(narnea.dist)
  narnea.clust.mrs <- kw_cluster_mrs(narnea.obj, narnea.clust$opt.clust)
  saveRDS(narnea.dist, file = paste('anchor-combined_analysis/', ct.name, '/',
                                    ct.name, '_int-seurat_vs-vehicle-ges_narnea-dist.rds', sep = ''))
  saveRDS(narnea.mds, file = paste('anchor-combined_analysis/', ct.name, '/',
                                   ct.name, '_int-seurat_vs-vehicle-ges_narnea-mds.rds', sep = ''))
  saveRDS(narnea.clust, file = paste('anchor-combined_analysis/', ct.name, '/',
                                     ct.name, '_int-seurat_vs-vehicle-ges_narnea-clust.rds', sep = ''))
  saveRDS(narnea.clust.mrs, file = paste('anchor-combined_analysis/', ct.name, '/',
                                         ct.name, '_int-seurat_vs-vehicle-ges_narnea-clust-mrs.rds', sep = ''))
}
###############

#############################################


### group vehicle networks
###############
sample.names <- c('KA003', 'KA005', 'KA006', 'KA007')
filt.counts <- readRDS('combined_analysis/combined-gexp_filt-counts.rds')
seurat.obj <- readRDS('combined_analysis/combined-gexp_seurat-obj.rds')
cluster.mappings <- read.csv('combined_analysis/cluster-mapping_20221107.csv')
colnames(cluster.mappings) <- c('Cluster', 'Cell.Type', 'Sub.Type')
## create label vectors
sample.name.vec <- sapply(colnames(seurat.obj), function(x) {strsplit(x, '\\.')[[1]][1]} )
group.name.vec <- sample.table$Group[match(sample.name.vec, sample.table$Sample)]
names(group.name.vec) <- names(sample.name.vec)
clust.vec <- seurat.obj$seurat_clusters
ct.vec <- plyr::mapvalues(clust.vec, from = cluster.mappings$Cluster, to = cluster.mappings$Cell.Type)
## table of vehicle samps and cell types
vehicle.samps <- names(group.name.vec)[which(group.name.vec == 'Vehicle')]
vehicle.ct.table <- table(sample.name.vec[vehicle.samps], ct.vec[vehicle.samps])
## make aracne matrices
dir.create('combined_analysis/a3_mats/')
for (gn in unique(ct.vec)) {
  print(gn)
  # make matrix
  gn.samps <- intersect(vehicle.samps, names(ct.vec)[which(ct.vec == gn)])
  a3.mat <- filt.counts[,gn.samps]
  a3.mat <- a3.mat[which(rowSums(a3.mat) > 0),]
  a3.mat <- cpm_norm(a3.mat, l2 = TRUE)
  if (ncol(a3.mat) > 1000) {a3.mat <- a3.mat[, sample(gn.samps, 1000)]}
  print(dim(a3.mat))
  # save
  gn.name <- tolower(gn)
  saveRDS(a3.mat, paste('combined_analysis/a3_mats/', gn.name, '_a3-mat.rds', sep = ''))
}
###############

### create network lists
###############
network.names <- c('endo', 'epi', 'fibro', 'lymphocyte', 'myeloid')
dir.create('combined_analysis/pad-a3_nets/net-lists/')
network.mdata <- list()
for (nn in network.names) {
  print(nn)
  nn.tf <- readRDS(paste('combined_analysis/pad-a3_nets/ad-', nn, '_tf_regulon-list.rds', sep = ''))
  nn.sig <- readRDS(paste('combined_analysis/pad-a3_nets/ad-', nn, '_sig_regulon-list.rds', sep = ''))
  nn.surf <- readRDS(paste('combined_analysis/pad-a3_nets/ad-', nn, '_surf_regulon-list.rds', sep = ''))
  nn.net <- c(nn.tf, nn.sig, nn.surf)
  saveRDS(nn.net, file = paste('combined_analysis/pad-a3_nets/net-lists/ad-', nn, '_net-list.rds', sep = ''))
  mdata.vec <- c(length(nn.tf), length(nn.sig), length(nn.surf), length(nn.net))
  network.mdata[[nn]] <- mdata.vec
}
mdata.df <- as.data.frame(Reduce(rbind, network.mdata))
rownames(mdata.df) <- network.names
colnames(mdata.df) <- c('TF', 'Sig', 'Surf', 'Total')
write.csv(mdata.df, file = 'combined_analysis/pad-a3_nets/net-size_metadata.csv', quote = FALSE)
###############

### make gene-union combined count matrix
###############
sample.names <- c('KA003', 'KA005', 'KA006', 'KA007')
sample.inds <- match(sample.names, sample.table$Sample)
## load count matrices
sample.counts <- list()
for (sn in sample.names) {
  print(sn)
  c.mat <- readRDS(paste('qc_expression/', sn, '_filt-counts.rds', sep = ''))
  sample.counts[[sn]] <- c.mat
}
## get union genes
union.genes <- unique(unlist(sapply(sample.counts, rownames)))
## pad each matrix with missing genes
pad.sample.counts <- lapply(sample.counts, function(x) {
  missing.genes <- setdiff(union.genes, rownames(x))
  pad.mat <- matrix(0L, nrow = length(missing.genes), ncol = ncol(x))
  rownames(pad.mat) <- missing.genes
  x.mat <- rbind(x, pad.mat)
  x.mat <- x.mat[union.genes,]
  return(x.mat)
})
## modify column names
sample.name.pref <- rep(sample.names, sapply(pad.sample.counts, ncol))
## combine into final matrix
pad.counts.mat <- Reduce(cbind, pad.sample.counts)
colnames(pad.counts.mat) <- paste(sample.name.pref, colnames(pad.counts.mat), sep = '.')
saveRDS(pad.counts.mat, file = 'combined_analysis/combined-gexp_pad-filt-counts.rds')
###############

### group vehicle networks w/ padded counts
###############
sample.names <- c('KA003', 'KA005', 'KA006', 'KA007')
filt.counts <- readRDS('combined_analysis/combined-gexp_pad-filt-counts.rds')
seurat.obj <- readRDS('combined_analysis/combined-gexp_seurat-obj.rds')
cluster.mappings <- read.csv('combined_analysis/cluster-mapping_20221107.csv')
colnames(cluster.mappings) <- c('Cluster', 'Cell.Type', 'Sub.Type')
## create label vectors
sample.name.vec <- sapply(colnames(seurat.obj), function(x) {strsplit(x, '\\.')[[1]][1]} )
group.name.vec <- sample.table$Group[match(sample.name.vec, sample.table$Sample)]
names(group.name.vec) <- names(sample.name.vec)
clust.vec <- seurat.obj$seurat_clusters
ct.vec <- plyr::mapvalues(clust.vec, from = cluster.mappings$Cluster, to = cluster.mappings$Cell.Type)
## table of vehicle samps and cell types
vehicle.samps <- names(group.name.vec)[which(group.name.vec == 'Vehicle')]
vehicle.ct.table <- table(group.name.vec[vehicle.samps], ct.vec[vehicle.samps])
## make aracne matrices
dir.create('combined_analysis/pad-a3_mats/')
for (gn in unique(ct.vec)) {
  print(gn)
  # make matrix
  gn.samps <- intersect(vehicle.samps, names(ct.vec)[which(ct.vec == gn)])
  a3.mat <- filt.counts[,gn.samps]
  a3.mat <- a3.mat[which(rowSums(a3.mat) > 0),]
  a3.mat <- cpm_norm(a3.mat, l2 = TRUE)
  if (ncol(a3.mat) > 1000) {a3.mat <- a3.mat[, sample(gn.samps, 1000)]}
  print(dim(a3.mat))
  # save
  gn.name <- tolower(gn)
  saveRDS(a3.mat, paste('combined_analysis/pad-a3_mats/', gn.name, '_a3-mat.rds', sep = ''))
}
###############

### group specific NaRnEA w/ Vehicle as reference
###############
counts.mat <- readRDS('combined_analysis/combined-gexp_filt-counts.rds')
seurat.obj <- readRDS('combined_analysis/combined-gexp_seurat-obj.rds')
cluster.mappings <- read.csv('combined_analysis/cluster-mapping_20221107.csv')
colnames(cluster.mappings) <- c('Cluster', 'Cell.Type', 'Sub.Type')
sample.name.vec <- sapply(colnames(seurat.obj), function(x) {strsplit(x, '\\.')[[1]][1]} )
group.name.vec <- sample.table$Group[match(sample.name.vec, sample.table$Sample)]
names(group.name.vec) <- names(sample.name.vec)
clust.vec <- seurat.obj$seurat_clusters
ct.vec <- plyr::mapvalues(clust.vec, from = cluster.mappings$Cluster, to = cluster.mappings$Cell.Type)
## map from group names to network names
name.map.df <- data.frame('CT' = unique(ct.vec),
                          'Network' = c('lymphocyte', 'myeloid', 'endo', 'epi', 'fibro'))
## PISCES analysis for each cell type group
for (ct in c('Lymphocyte', 'Endothelial', 'Epithelial', 'Fibroblast')) {
  print(ct)
  # set name and make directory
  ct.file.name <- tolower(ct)
  dir.create(paste('combined_analysis/narnea_analysis/', ct.file.name, sep = ''))
  
  # load network
  net.name <- name.map.df$Network[match(ct, name.map.df$CT)]
  net.list <- readRDS(paste('combined_analysis/a3_nets/net-lists/ad-', net.name, '_net-list.rds', sep = ''))
  
  # get cells for this cell type
  ct.cells <- names(ct.vec)[which(ct.vec == ct)]
  ct.counts <- counts.mat[, ct.cells]
  ct.counts <- ct.counts[which(rowSums(ct.counts) > 0),]
  ct.counts <- cpm_norm(ct.counts, l2 = TRUE)
  saveRDS(ct.counts, file = paste('combined_analysis/narnea_analysis/', ct.file.name, '/', 
                                  ct.file.name, '_cpm-counts.rds', sep = ''))
  
  # generate signature against vehicle
  vehicle.cells <- names(which(group.name.vec[ct.cells] == 'Vehicle'))
  vehicle.mean <- apply(ct.counts[, vehicle.cells], 1, mean)
  vehicle.sd <- apply(ct.counts[, vehicle.cells], 1, sd)
  use.genes <- which(vehicle.sd != 0)
  ct.ges <- apply(ct.counts[use.genes,], 2, function(x) {
    (x - vehicle.mean[use.genes]) / vehicle.sd[use.genes]
  })
  saveRDS(ct.ges, file = paste('combined_analysis/narnea_analysis/', ct.file.name, '/',
                               ct.file.name, '_vs-vehicle-ges.rds', sep = ''))
  
  # run narnea
  ct.narnea <- matrix_narnea(ct.ges, net.list)
  saveRDS(ct.narnea, file = paste('combined_analysis/narnea_analysis/', ct.file.name, '/',
                                  ct.file.name, '_vs-vehicle-narnea.rds', sep = ''))
  
  # narnea objects
  narnea.dist <- dist(t(ct.narnea$PES))
  narnea.mds <- make_mds(narnea.dist)
  narnea.clust <- louvain_k(narnea.dist)
  narnea.clust.mrs <- kw_cluster_mrs(ct.narnea, narnea.clust$opt.clust)
  saveRDS(narnea.mds, file = paste('combined_analysis/narnea_analysis/', ct.file.name, '/',
                                   ct.file.name, '_vs-vehicle-narnea_mds.rds', sep = ''))
  saveRDS(narnea.dist, file = paste('combined_analysis/narnea_analysis/', ct.file.name, '/',
                                    ct.file.name, '_vs-vehicle-narnea_dist.rds', sep = ''))
  saveRDS(narnea.clust, file = paste('combined_analysis/narnea_analysis/', ct.file.name, '/',
                                     ct.file.name, '_vs-vehicle-narnea_clust.rds', sep = ''))
  saveRDS(narnea.clust.mrs, file = paste('combined_analysis/narnea_analysis/', ct.file.name, '/',
                                         ct.file.name, '_vs-vehicle-narnea_clust-mrs.rds', sep = ''))
}
## make plots for each cell type
for (ct in c('Lymphocyte', 'Endothelial', 'Epithelial', 'Fibroblast')) {
  print(ct)
  ct.file.name <- tolower(ct)
  
  # load objects
  ct.narnea <- readRDS(paste('combined_analysis/narnea_analysis/', ct.file.name, '/',
                             ct.file.name, '_vs-vehicle-narnea.rds', sep = ''))
  narnea.mds <- readRDS(paste('combined_analysis/narnea_analysis/', ct.file.name, '/',
                              ct.file.name, '_vs-vehicle-narnea_mds.rds', sep = ''))
  narnea.clust <- readRDS(paste('combined_analysis/narnea_analysis/', ct.file.name, '/',
                                ct.file.name, '_vs-vehicle-narnea_clust.rds', sep = ''))
  narnea.clust.mrs <- readRDS(paste('combined_analysis/narnea_analysis/', ct.file.name, '/',
                                    ct.file.name, '_vs-vehicle-narnea_clust-mrs.rds', sep = ''))
  ct.cells <- names(ct.vec)[which(ct.vec == ct)]
  
  ## scatter plots
  plot.df <- data.frame('MDS1' = narnea.mds[ct.cells, 1], 'MDS2' = narnea.mds[ct.cells, 2],
                        'Sample' = sample.name.vec[ct.cells],
                        'Group' = group.name.vec[ct.cells],
                        'Cluster' = as.factor(narnea.clust$opt.clust[ct.cells]))
  # calculate centers
  cluster.mds.centers <- lapply(sort(unique(narnea.clust$opt.clust)), function(x) {
    clust.samps <- which(narnea.clust$opt.clust == x)
    mean.vals <- colMeans(narnea.mds[clust.samps,])
    return(mean.vals)
  })
  cluster.mds.centers <- Reduce(rbind, cluster.mds.centers)
  cluster.mds.centers <- as.data.frame(cluster.mds.centers)
  cluster.mds.centers$Cluster <- sort(unique(narnea.clust$opt.clust))
  # sample
  ggplot(plot.df, aes(MDS1, MDS2)) + geom_point(aes(color = Sample)) + 
    ggtitle(paste(ct, ' Activity Clustering - Sample', sep = '')) + 
    guides(colour = guide_legend(override.aes = list(size = 3)))
  ggsave(paste('combined_analysis/narnea_analysis/', ct.file.name, '/',
               ct.file.name, '_vs-vehicle-narnea_sample-mds.jpg', sep = ''),
         height = 8, width = 10, units = 'in', dpi = 300)
  # group
  ggplot(plot.df, aes(MDS1, MDS2)) + geom_point(aes(color = Group)) + 
    ggtitle(paste(ct, ' Activity Clustering - Group', sep = '')) + 
    guides(colour = guide_legend(override.aes = list(size = 3)))
  ggsave(paste('combined_analysis/narnea_analysis/', ct.file.name, '/',
               ct.file.name, '_vs-vehicle-narnea_group-mds.jpg', sep = ''),
         height = 8, width = 10, units = 'in', dpi = 300)
  # cluster
  ggplot(plot.df, aes(MDS1, MDS2)) + geom_point(aes(color = Cluster)) + 
    ggtitle(paste(ct, ' Activity Clustering - Cluster', sep = '')) + 
    guides(colour = guide_legend(override.aes = list(size = 3))) + 
    annotate("text", x = cluster.mds.centers$V1, y = cluster.mds.centers$V2,
             label = cluster.mds.centers$Cluster)
  ggsave(paste('combined_analysis/narnea_analysis/', ct.file.name, '/',
               ct.file.name, '_vs-vehicle-narnea_cluster-mds.jpg', sep = ''),
         height = 8, width = 10, units = 'in', dpi = 300)
  
  # regulator heamtap
  ct.heatmap <- cluster_mr_heatmap(ct.narnea$PES, dat.type = 'pact', clust.vec = narnea.clust$opt.clust,
                                   mr.list = narnea.clust.mrs, reg.class = 'regulator', scale.rows = TRUE,
                                   plot.title = paste(ct, ' Cluster MRs (Regulators)', sep = ''))
  jpeg(paste('combined_analysis/narnea_analysis/', ct.file.name, '/',
             ct.file.name, '_vs-vehicle-narnea_cluster-mrs_regs.jpg', sep = ''),
       height = 10, width = 8, res = 300, units = 'in')
  print(ct.heatmap)
  dev.off()
  # marker heatmap
  ct.heatmap <- cluster_mr_heatmap(ct.narnea$PES, dat.type = 'pact', clust.vec = narnea.clust$opt.clust,
                                   mr.list = narnea.clust.mrs, reg.class = 'marker', scale.rows = TRUE,
                                   plot.title = paste(ct, ' Cluster MRs (Markers)', sep = ''))
  jpeg(paste('combined_analysis/narnea_analysis/', ct.file.name, '/',
             ct.file.name, '_vs-vehicle-narnea_cluster-mrs_markers.jpg', sep = ''),
       height = 10, width = 8, res = 300, units = 'in')
  print(ct.heatmap)
  dev.off()
}
###############

### endothelial tip marker plots
###############
tip.markers <- read.csv('pdac_endothelial-markers.csv')
colnames(tip.markers) <- c('Gene', 'Group')
tip.markers <- tip.markers[order(tip.markers$Group),]
## load narnea objects
narnea.res <- readRDS('combined_analysis/narnea_analysis/endothelial/endothelial_vs-vehicle-narnea.rds')
narnea.clust <- readRDS('combined_analysis/narnea_analysis/endothelial/endothelial_vs-vehicle-narnea_clust.rds')
clust.vec <- narnea.clust$opt.clust
## pull out PES mat and remove NAs (TO DO: check why there are NAs)
pes.mat <- narnea.res$PES
row.has.na <- apply(pes.mat, 1, function(x) {length(which(is.na(x))) > 0})
pes.mat <- pes.mat[-which(row.has.na),]
## create annotation vectors
sample.name.vec <- sapply(colnames(pes.mat), function(x) {strsplit(x, '\\.')[[1]][1]} )
group.name.vec <- sample.table$Group[match(sample.name.vec, sample.table$Sample)]
names(group.name.vec) <- names(sample.name.vec)
## make plot
cell.order <- names(sort(clust.vec))
plot.mat <- pes.mat[, cell.order]
marker.set <- unique(intersect(rownames(plot.mat), tip.markers$Gene))
marker.group <- tip.markers$Group[match(marker.set, tip.markers$Gene)]
plot.mat <- plot.mat[marker.set,]
# set colors
col.breaks <- quantile_breaks(plot.mat)
col.fun <- circlize::colorRamp2(col.breaks, color_levels('pact', length(col.breaks)))
clust.colors <- group_colors(length(unique(clust.vec)))
names(clust.colors) <- sort(unique(clust.vec))
sample.colors <- group_colors(length(unique(sample.name.vec)))
names(sample.colors) <- sort(unique(sample.name.vec))
group.colors <- c('darkgrey', 'lightgrey')
names(group.colors) <- sort(unique(group.name.vec))
marker.group.colors <- c('green', 'orange')
names(marker.group.colors) <- sort(unique(tip.markers$Group))
# column annotation
column.annot <- HeatmapAnnotation('Cluster' = clust.vec[cell.order],
                                  'Sample' = sample.name.vec[cell.order],
                                  'Group' = group.name.vec[cell.order],
                                  col = list('Cluster' = clust.colors,
                                             'Sample' = sample.colors, 'Group' = group.colors))
col.gaps <- clust.vec[cell.order]
# row annotation
row.annot <- rowAnnotation('SubGroup' = marker.group,
                           col = list('SubGroup' = marker.group.colors), 
                           show_annotation_name = FALSE,
                           show_legend = TRUE)
row.gaps <- marker.group
# make plot
jpeg('combined_analysis/narnea_analysis/endothelial/endothelial_vs-vehicle-narnea_stalk-tip-markers.jpg',
     height = 10, width = 15, units = 'in', res = 300)
print(Heatmap(plot.mat, name = 'PES',
              col = col.fun,
              top_annotation = column.annot, column_split = col.gaps,
              left_annotation = row.annot, row_split = row.gaps,
              cluster_rows = FALSE, cluster_columns = FALSE,
              show_row_names = TRUE, show_column_names = FALSE,
              column_title = 'Endothelial Cells - Stalk/Tip Markers', row_title = NULL,
              row_names_gp = gpar(fontsize = 8)))
dev.off()
###############

### group-level fibroblast analysis
###############
counts.mat <- readRDS('combined_analysis/combined-gexp_pad-filt-counts.rds')
seurat.obj <- readRDS('combined_analysis/combined-gexp_seurat-obj.rds')
cluster.mappings <- read.csv('combined_analysis/cluster-mapping_20221107.csv')
colnames(cluster.mappings) <- c('Cluster', 'Cell.Type', 'Sub.Type')
sample.name.vec <- sapply(colnames(seurat.obj), function(x) {strsplit(x, '\\.')[[1]][1]} )
group.name.vec <- sample.table$Group[match(sample.name.vec, sample.table$Sample)]
names(group.name.vec) <- names(sample.name.vec)
clust.vec <- seurat.obj$seurat_clusters
ct.vec <- plyr::mapvalues(clust.vec, from = cluster.mappings$Cluster, to = cluster.mappings$Cell.Type)
## process matrices
cpm.mat <- cpm_norm(counts.mat, l2 = TRUE)
sct.mat <- as.matrix(seurat.obj@assays$SCT@scale.data)
## make expression boxplots
plot.df <- data.frame('Cell.Type' = ct.vec,
                      'Group' = group.name.vec)
plot.df <- cbind(plot.df, t(cpm.mat[fibro.genes,]))
plot.df <- melt(plot.df, id.vars = c('Cell.Type', 'Group'))
colnames(plot.df) <- c('Cell.Type', 'Group', 'Variable', 'L2CPM')
ggplot(plot.df, aes(x = Cell.Type, y = L2CPM)) + geom_boxplot(aes(color = Cell.Type)) +
  facet_grid(rows = vars(Variable), cols = vars(Group)) +
  ggtitle('Expression of Fibroblast Genes by Group')
ggsave('combined_analysis/narnea_analysis/fibroblast/fibroblast_exp-boxplot.jpg',
       height = 6, width = 6, units = 'in', dpi = 300)
## get ct groups
ct <- 'Fibroblast'
ct.net <- readRDS('combined_analysis/a3_nets/net-lists/ad-fibro_net-list.rds')
ct.cells <- names(ct.vec)[which(ct.vec == 'Fibroblast')]
vehicle.samps <- names(group.name.vec)[which(group.name.vec == 'Vehicle')]
treatment.samps <- names(group.name.vec)[which(group.name.vec == '926')]
ct.vehicle.cells <- intersect(ct.cells, vehicle.samps)
ct.926.cells <- intersect(ct.cells, treatment.samps)
## make signature of 926 vs vehicle
ges.vec <- apply(cpm.mat, 1, function(x) {
  w.test <- wilcox.test(x[ct.926.cells], x[ct.vehicle.cells], alternative = "two.sided")
  p.val <- w.test$p.val
  rbs.cor <- 2 * w.test$statistic / (length(ct.926.cells) * length(ct.vehicle.cells)) - 1
  ges.val <- qnorm(p.val, lower.tail = FALSE) * sign(rbs.cor)
  return(ges.val)
})
## run narnea
if(length(which(is.na(ges.vec))) > 0) {ges.vec <- ges.vec[-which(is.na(ges.vec))]}
ges.mat <- matrix(ges.vec, nrow = length(ges.vec))
rownames(ges.mat) <- names(ges.vec)
colnames(ges.mat) <- '926.v.vehicle'
narnea.res <- matrix_narnea(ges.mat, ct.net)
saveRDS(ges.vec, file = 'combined_analysis/narnea_analysis/fibroblast/fibroblast_926-v-vehicle_ges.rds')
saveRDS(narnea.res, file = 'combined_analysis/narnea_analysis/fibroblast/fibroblast_926-v-vehicle_narnea.rds')
###############

### fibroblast marker plots
###############
cpm.mat <- readRDS('combined_analysis/narnea_analysis/fibroblast/fibroblast_cpm-counts.rds')
narnea.res <- readRDS('combined_analysis/narnea_analysis/fibroblast/fibroblast_vs-vehicle-narnea.rds')
narnea.clust <- readRDS('combined_analysis/narnea_analysis/fibroblast/fibroblast_vs-vehicle-narnea_clust.rds')
clust.vec <- narnea.clust$opt.clust
## general heatmap variables
cell.order <- names(sort(clust.vec))
clust.colors <- group_colors(length(unique(clust.vec)))
names(clust.colors) <- sort(unique(clust.vec))
class.colors <- group_colors(length(unique(fibro.class.markers$Group)))
names(class.colors) <- unique(fibro.class.markers$Group)
column.annot <- HeatmapAnnotation('Cluster' = clust.vec[cell.order],
                                  col = list('Cluster' = clust.colors))
col.gaps <- clust.vec[cell.order]
## make heatmap with CAF marker activity
marker.use.inds <- which(fibro.class.markers$Gene %in% rownames(narnea.res$PES))
plot.mat <- narnea.res$PES[fibro.class.markers$Gene[marker.use.inds], cell.order]
plot.mat <- t(apply(plot.mat, 1, scale))
colnames(plot.mat) <- cell.order
col.breaks <- quantile_breaks(plot.mat)
col.fun <- circlize::colorRamp2(col.breaks, color_levels('pact', length(col.breaks)))
group.vec <- fibro.class.markers$Group[marker.use.inds]
names(group.vec) <- fibro.class.markers$Gene[marker.use.inds]
row.annot <- rowAnnotation('Group' = group.vec,
                           col = list('Group' = class.colors))
row.gaps <- group.vec
Heatmap(plot.mat, name = 'PES',
        col = col.fun,
        top_annotation = column.annot, left_annotation = row.annot,
        column_split = col.gaps, row_split = row.gaps,
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = FALSE,
        column_title = 'Fibroblast PISCES Clustering - Marker Activity', row_title = NULL,
        row_names_gp = gpar(fontsize = 8))
## make heatmap with CAF marker expression
marker.use.inds <- which(fibro.class.markers$Gene %in% rownames(cpm.mat))
plot.mat <- cpm.mat[fibro.class.markers$Gene[marker.use.inds], cell.order]
plot.mat <- t(apply(plot.mat, 1, scale))
colnames(plot.mat) <- cell.order
col.breaks <- quantile_breaks(plot.mat)
col.fun <- circlize::colorRamp2(col.breaks, color_levels('gexp', length(col.breaks)))
group.vec <- fibro.class.markers$Group[marker.use.inds]
names(group.vec) <- fibro.class.markers$Gene[marker.use.inds]
row.annot <- rowAnnotation('Group' = group.vec,
                           col = list('Group' = class.colors))
row.gaps <- group.vec
Heatmap(plot.mat, name = 'L2CPM',
        col = col.fun,
        top_annotation = column.annot, left_annotation = row.annot,
        column_split = col.gaps, row_split = row.gaps,
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = FALSE,
        column_title = 'Fibroblast PISCES Clustering - Marker Expression', row_title = NULL,
        row_names_gp = gpar(fontsize = 8))
###############



