# Examining transcriptional waves using RNA velocity through IDC

# Colours
require("RColorBrewer")

## 10X
dat_10X <- readRDS("bergei_10X_velodata.rds")
velo_10X <- readRDS("bergei_10X_veloEst.rds")
pca_10X <- readRDS("bergei_10X_veloPCA.rds")


# Set up SCE
require("scater")
require("scmap")
sce <- SingleCellExperiment(assays=list(counts=velo_10X$current, projected_counts=velo_10X$projected), colData=dat_10X$ann)
rowData(sce)$feature_symbol <- rownames(assays(sce)[["counts"]])
sce <- calculateQCMetrics(sce)
sce <- normalise(sce)
assays(sce)[["logcounts"]] <- as.matrix(assays(sce)[["logcounts"]])
assays(sce)[["counts"]] <- as.matrix(assays(sce)[["counts"]])


proj_mat <- assays(sce)[["logcounts"]] + 5*velo_10X$deltaE
proj_mat[proj_mat < 0] <- 0
sce_p <- SingleCellExperiment(assays=list(logcounts=as.matrix(proj_mat)), colData=dat_10X$ann)
rowData(sce_p)$feature_symbol <- rownames(assays(sce_p)[["logcounts"]])

# Build index
sce <- selectFeatures(sce, suppress_plot = FALSE, n_features = 500)
set.seed(1)
sce <- indexCell(sce)

# Map projected ### THIS DOESN'T REALLY WORK NEED OTHER PCA MATRIX
scmapCell_results <- scmapCell(
  sce_p, 
  list(
    current = metadata(sce)$scmap_cell_index
  ),
  w=10
)


get_pcs <- function(a) {
	Matrix::colMeans(as.matrix(colData(sce)[a,c("pc1", "pc2")]))
}

projected_pcs <- t(apply(scmapCell_results$current$cells, 2, get_pcs))
current_pcs <- as.matrix(colData(sce)[,c("pc1", "pc2")])


plot(current_pcs, pch=16, col=dat_10X$ann$cell_color)
arrows(current_pcs[,1], current_pcs[,2], projected_pcs[,1], projected_pcs[,2], len=0.1, lty=3)


gridn=30
gridxes <- seq(from=min(projected_pcs[,1], current_pcs[,1])-0.0001, to=max(projected_pcs[,1], current_pcs[,1])+0.0001, length=gridn)
gridyes <- seq(from=min(projected_pcs[,2], current_pcs[,2])-0.0001, to=max(projected_pcs[,2], current_pcs[,2])+0.0001, length=gridn)

gridpts <- cbind(rep(gridxes, times=gridn), rep(gridyes, each=gridn))
binned <- cbind(cut(current_pcs[,1], gridxes), cut(current_pcs[,2], gridyes))

par(mar=c(4,4,1,1))
plot(current_pcs, pch=21, bg=dat_10X$ann$cell_color, col="grey65", cex=1)
points(rep(gridxes, times=gridn), rep(gridyes, each=gridn), pch=16, cex=0.3)

for(i in min(binned[,1]):max(binned[,1])) {
for(j in min(binned[,1]):max(binned[,1])) {
        thispts <- binned[,1] == i & binned[,2] == j
        if (sum(thispts) < 5) {next;}
        arrows( mean(current_pcs[thispts,1]),
                mean(current_pcs[thispts,2]),
                mean(projected_pcs[thispts,1]),
                mean(projected_pcs[thispts,2]), len=0.07 )

}
}

