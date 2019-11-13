# Read in existing PCA and cell labels:
pca1 <- read.delim("/nfs/users/nfs_t/ta6/Collaborations/MCA/For10XVelocity/mca_pb_tmm_PCA.csv", ",", header=TRUE) #columns: cellID, pc1, pc2
labels <- read.delim("/nfs/users/nfs_t/ta6/Collaborations/MCA/For10XVelocity/stage_pred.csv", ",", header=TRUE)
labels <- labels[match(pca1$cell_id, labels[,1]),]
labels <- labels[!is.na(labels[,1]),]

# Read in velocity counts
emat <- readRDS("out_S_Matrix.rds")
smat <- readRDS("out_U_Matrix.rds")

require("Matrix")

# read in and clean up cell names
cell_names <- read.delim("out_ca_matrix.csv", ",", header=TRUE)
cell_names <- sub("possorted_genome_bam_UVJ4K:", "", cell_names$CellID)
cell_names <- sub("x$", "", cell_names)

# read in gene names
gene_ann <- read.delim("out_ra_matrix.csv", ",", header=TRUE)
gene_names <- gene_ann$Accession

rownames(emat) <- gene_names
rownames(smat) <- gene_names
colnames(emat) <- cell_names
colnames(smat) <- cell_names

# Cell and Gene QC

# Consistent cells
emat <- emat[,match(labels[,1], colnames(emat))]
smat <- smat[,match(labels[,1], colnames(smat))]

# Expressed genes
keep_g <- rowSums(emat > 0) > 20 & rowSums(smat > 0) > 10

emat <- emat[keep_g,] 
smat <- smat[keep_g,] 

# combine annotations and set up known cell-type colours
source("/nfs/users/nfs_t/ta6/Collaborations/MCA/For10XVelocity/Colour_Scheme.R")
ann <- labels
rownames(ann) <- labels[,1]
pca1 <- pca1[match(labels[,1], pca1$cell_id),]
ann$pc1 <- pca1[,2]
ann$pc2 <- pca1[,3]
ann$cell_id <- pca1[,1]

colors <- colors[order(as.numeric(names(colors)))]
names <- names[order(as.numeric(names(names)))]
ann$cell_color <- colors[labels[,2]+1]
ann$cell_name <- names[labels[,2]+1]

# Save Nice input
dat <- list(emat=emat, smat=smat, ann=ann)
saveRDS(dat, file="bergei_10X_velodata.rds")

### Fitting ###

dat_local <- dat
require("velocyto.R")

set.seed(4724)
rvel.qf <- gene.relative.velocity.estimates(dat_local$emat, dat_local$smat, deltaT=1, kCells = 50, fit.quantile = 0.2, min.nmat.emat.slope=0.2, min.nmat.emat.correlation=0.2, kGenes=5)

# Exclude poorly fit genes

poor_fits <- apply(rvel.qf$deltaE, 1, function(x) {x <- x[x!=0]; a <- table(sign(x)); max(a/sum(a))})
poor_fit_genes <- names(poor_fits)[poor_fits > 0.90]

dat_local$emat <- dat_local$emat[! rownames(dat_local$emat) %in% poor_fit_genes,]
dat_local$smat <- dat_local$smat[! rownames(dat_local$smat) %in% poor_fit_genes,]

# rerun fitting

set.seed(4724)
rvel.qf <- gene.relative.velocity.estimates(dat_local$emat, dat_local$smat, deltaT=1, kCells = 50, fit.quantile = 0.2, min.nmat.emat.slope=0.2, min.nmat.emat.correlation=0.2, kGenes=5)


saveRDS(rvel.qf, file="bergei_10X_veloEst.rds")

### Velocity PCA ###
my_colors <- dat_local$ann$cell_color
names(my_colors) <- rownames(dat_local$ann)
png("bergei_10X_veloPCA.png", width=7, height=7, units="in", res=300)
pca.out <- pca.velocity.plot(rvel.qf,nPcs=5,plot.cols=2,cex=1, show.grid.flow=TRUE, grid.n=25, cell.colors=my_colors)
dev.off()
saveRDS(pca.out, file="bergei_10X_veloPCA.rds")


# Set-up transcription rate colours
require("RColorBrewer")
require("Matrix")
up_col <- brewer.pal(8, "Reds")
down_col <- brewer.pal(8, "Blues")

## 10X
dat_10X <- readRDS("bergei_10X_velodata.rds")
velo_10X <- readRDS("bergei_10X_veloEst.rds")

### Cell-level scores
# total transcription
tmp<- velo_10X$deltaE
tmp[tmp < 0] <- 0
cell_up <- Matrix::colSums(tmp)
# total decreasing expression
tmp<- velo_10X$deltaE
tmp[tmp > 0] <- 0
cell_down <- Matrix::colSums(tmp)
# net change
tmp<- velo_10X$deltaE
cell_both <- Matrix::colSums(tmp)

# Bin scores
breaks_down <- seq(from=min(cell_down), to=0, length=7)
breaks_up <- seq(from=0, to=max(cell_up), length=7)

cell_down_col <- rev(down_col)[cut(cell_down, breaks=breaks_down)]
cell_up_col <- up_col[cut(cell_up, breaks=breaks_up)]


### Plot ###
png("10X_IDC_Waves_Down.png", width=6, height=6, units="in", res=300)
plot(dat_10X$ann$pc1, dat_10X$ann$pc2, col=cell_down_col, pch=16, xlab="PC1", ylab="PC2")
dev.off()
png("10X_IDC_Waves_Up.png", width=6, height=6, units="in", res=300)
plot(dat_10X$ann$pc1, dat_10X$ann$pc2, col=cell_up_col, pch=16, xlab="PC1", ylab="PC2")
dev.off()
png("10X_IDC_Waves_Type.png", width=6, height=6, units="in", res=300)
plot(dat_10X$ann$pc1, dat_10X$ann$pc2, col=dat_10X$ann$cell_color, pch=16, xlab="PC1", ylab="PC2")
dev.off()

