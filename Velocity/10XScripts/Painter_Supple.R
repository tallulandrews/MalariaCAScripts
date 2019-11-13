# Read in published
extern <- read.delim("/nfs/users/nfs_t/ta6/Collaborations/MCA/Pf_transcription_vs_time_FromSuppl1.csv", sep=",", header=T)
extern <- extern[!is.na(extern[,51]), ]
stable <- read.delim("/nfs/users/nfs_t/ta6/Collaborations/MCA/Pf_stabilitzation_vs_time_FromSuppl1.csv", sep=",", header=T, stringsAsFactors=FALSE)
stable[,51] <- as.numeric(stable[,51])
stable <- stable[!is.na(stable[,51]), ]
ortho <- read.table("/nfs/users/nfs_t/ta6/Collaborations/MCA/orthologueTable.csv", sep=",", header=T)
ortho <- ortho[ortho[,4] %in% extern$Gene.ID,]
ortho[,6] <- sub("PKH", "PKNH", ortho[,6])
ortho[,6] <- paste(ortho[,6], "0", sep="")

## 10X velocyto
pb_dat_10X <- readRDS("bergei_10X_velodata_wPseudo.rds")
pb_velo_10X <- readRDS("bergei_10X_veloEst.rds")
pk_dat_10X <- readRDS("knowlesi_10X_velodata.rds")
pk_velo_10X <- readRDS("knowlesi_10X_veloEst.rds")
pf_dat_10X <- readRDS("falciparum_10X_velodata.rds")
pf_velo_10X <- readRDS("falciparum_10X_veloEst.rds")

pseudotime <- pb_dat_10X$ann$pseudotime
cell_color <- pb_dat_10X$ann$cells_recolor
# bin cells by pseudotime
n_pbins=20
p_breaks <- seq(from=0, to=2*pi, length=n_pbins+1)
p_bins <- cut(pseudotime, breaks=n_pbins)
p_axis_ats <- seq(from=0, to=1, length=n_pbins)
p_axis_labs <- round( (p_breaks[1:n_pbins]+p_breaks[2:(n_pbins+1)])/2, digits=1 )

# Get Orthos
extern <- extern[extern$Gene.ID %in% ortho[,4],]
extern$pbanka <- ortho[match(extern$Gene.ID, ortho[,4]), 1]
extern$pfalcip <- ortho[match(extern$Gene.ID, ortho[,4]), 4]
extern$pknowl <- ortho[match(extern$Gene.ID, ortho[,4]), 6]
stable <- stable[stable$Gene.ID %in% ortho[,4],]
stable$pbanka <- ortho[match(stable$Gene.ID, ortho[,4]), 1]
stable$pfalcip <- ortho[match(stable$Gene.ID, ortho[,4]), 4]
stable$pknowl <- ortho[match(stable$Gene.ID, ortho[,4]), 6]

# Match external
up_mat<- pb_velo_10X$deltaE
up_mat[up_mat < 0] <- 0
cell_up <- Matrix::colSums(up_mat)

spliced_mat <- pb_dat_10X$emat
spliced_mat <- spliced_mat[Matrix::rowSums(spliced_mat>0)>20,]

spp <- "pbanka"

up_mat <- up_mat[rownames(up_mat) %in% extern[,spp], ]
spliced_mat <- spliced_mat[rownames(spliced_mat) %in% extern[,spp], ]
extern_spp <- extern[match(rownames(up_mat), extern[,spp]), ]
stable_spp <- stable[match(rownames(spliced_mat), stable[,spp]), ]

# Colours
require("RColorBrewer")
require("Matrix")
up_col <- brewer.pal(8, "Reds")
down_col <- brewer.pal(8, "Blues")

# aggregate
fixed_row_mean_aggregate <-function (mat, groups) {
    MAT <- as.matrix(mat)
    x <- split(seq(ncol(MAT)), groups)
    result <- sapply(x, function(a) {
                        if (length(a) > 1) {
                        rowMeans(MAT[, a])
                        } else if (length(a) == 1) {
                        MAT[,a]
                        } else {
                        rep(0, times=nrow(MAT))
                        }
                })
    return(result)
}

# bin genes by pf peak time
n_gbins=36
g_breaks <- seq(from=1, to=48, length=n_gbins+1)
tg_bins <- cut(extern_spp[,51], breaks=n_gbins);
sg_bins <- cut(stable_spp[,51], breaks=n_gbins);
g_axis_ats <- seq(from=0, to=1, length=n_gbins)
g_axis_labs <- signif( (g_breaks[1:n_gbins]+g_breaks[2:(n_gbins+1)])/2, digits=2 )

# Scale t-rate
up_s_mat <- t(apply(up_mat, 1, function(a) {a<- a-min(a); a <- a/max(a);}))
cell_up_s <- Matrix::colSums(up_s_mat)
heat_data_s <- fixed_row_mean_aggregate(t(up_s_mat), tg_bins);
heat_data_s <- fixed_row_mean_aggregate(t(heat_data_s), p_bins);
# Scale spliced
spliced_s_mat <- t(apply(spliced_mat, 1, function(a) {a<- a-min(a); a <- a/max(a);}))
cell_spliced_s <- Matrix::colSums(spliced_s_mat)
heat_data_e <- fixed_row_mean_aggregate(t(spliced_s_mat), sg_bins);
heat_data_e <- fixed_row_mean_aggregate(t(heat_data_e), p_bins);

# Supplementary Figure
png("bergei_SupplFig_Scaled_TranscriptionalWaves_10X.png", width=8*2, height=5, units="in", res=300)
graphics::layout(mat=rbind(c(8,5,1,4), c(7,6,2,3)), widths=c(2,5,5,2), heights=c(2,5))
par(mar=c(0,4,1,1))

# trans through time
plot(pseudotime, cell_up_s, xlab="pseudotime", ylab="increase in transcription", xlim=c(min(pseudotime)+0.22, max(pseudotime)-0.22), xaxt="n", col=cell_color)
abline(v=p_breaks, lty=3)
points(pseudotime, cell_up_s,col=cell_color)

# trans peaks Matching
par(mar=c(4,4,1,1))
image(x=p_breaks, y=g_breaks, z=t(heat_data_s), ylab="peak time (h)", xlab="pseudotime", col=up_col, yaxt="n")
axis(2, at=g_axis_labs[seq(from=1, to=length(g_axis_labs), by=5)], labels=g_axis_labs[seq(from=1, to=length(g_axis_labs), by=5)])

# trans gene Histogram
par(mar=c(4,0,1,1))
hist_vals <- table(tg_bins)
barplot(hist_vals, axes=FALSE, ylim=c(0+1.5, length(hist_vals)-1.5), xlim=c(0, max(hist_vals)), space=0, horiz=TRUE, border=NA, col="grey70", names="", xlab="nGenes", bty="l")
axis(side=1)
axis(side=2, at=c(-1,seq(from=1, to=length(g_axis_labs), by=5)-0.5,60), 
	labels=rep("", times=length(seq(from=1, to=length(g_axis_labs), by=5))+2))

# blank R corner
source("~/R-Scripts/Blank_plot.R")
blank_plot()

# expr through time
par(mar=c(0,1,1,4))
plot(pseudotime, cell_spliced_s, xlab="", ylab="", 
	xlim=c(min(pseudotime)+0.22, max(pseudotime)-0.22), 
	xaxt="n", yaxt="n", col=cell_color,
	ylim=c(0,150))
axis(side=4)
abline(v=p_breaks, lty=3)
points(pseudotime, cell_spliced_s,col=cell_color)

# spliced peaks Matching
par(mar=c(4,1,1,4))
image(x=p_breaks, y=g_breaks, z=t(heat_data_e), ylab="", xlab="pseudotime", col=down_col, yaxt="n")
axis(4, at=g_axis_labs[seq(from=1, to=length(g_axis_labs), by=5)], labels=g_axis_labs[seq(from=1, to=length(g_axis_labs), by=5)])

# trans gene Histogram
par(mar=c(4,1,1,0))
hist_vals <- table(sg_bins)
barplot(-hist_vals, axes=FALSE, ylim=c(0+1.5, length(hist_vals)-1.5), xlim=c(-max(hist_vals), 0), space=0, horiz=TRUE, border=NA, col="grey70", names="", xlab="nGenes", bty="l")
axis(side=1, at=seq(from=0, to=-100, by=-20), labels=seq(from=0, to=-100, by=-20)*-1)
axis(side=4, at=c(-1,seq(from=1, to=length(g_axis_labs), by=5)-0.5,60), 
	labels=rep("", times=length(seq(from=1, to=length(g_axis_labs), by=5))+2))

dev.off()

#### UNSCALED ####

# Scale t-rate
up_s_mat <- up_mat
cell_up_s <- Matrix::colSums(up_s_mat)
heat_data_s <- fixed_row_mean_aggregate(t(up_s_mat), tg_bins);
heat_data_s <- fixed_row_mean_aggregate(t(heat_data_s), p_bins);
# Scale spliced
spliced_s_mat <- spliced_mat
cell_spliced_s <- Matrix::colSums(spliced_s_mat)
heat_data_e <- fixed_row_mean_aggregate(t(spliced_s_mat), sg_bins);
heat_data_e <- fixed_row_mean_aggregate(t(heat_data_e), p_bins);

# Supplementary Figure
png("bergei_SupplFig_UnScaled_TranscriptionalWaves_10X.png", width=8*2, height=5, units="in", res=300)
graphics::layout(mat=rbind(c(8,5,1,4), c(7,6,2,3)), widths=c(2,5,5,2), heights=c(2,5))
par(mar=c(0,4,1,1))

# trans through time
plot(pseudotime, cell_up_s, xlab="pseudotime", ylab="increase in transcription", xlim=c(min(pseudotime)+0.22, max(pseudotime)-0.22), xaxt="n", col=cell_color)
abline(v=p_breaks, lty=3)
points(pseudotime, cell_up_s,col=cell_color)

# trans peaks Matching
par(mar=c(4,4,1,1))
image(x=p_breaks, y=g_breaks, z=t(log2(heat_data_s+1)), ylab="peak time (h)", xlab="pseudotime", col=up_col, yaxt="n")
axis(2, at=g_axis_labs[seq(from=1, to=length(g_axis_labs), by=5)], labels=g_axis_labs[seq(from=1, to=length(g_axis_labs), by=5)])

# trans gene Histogram
par(mar=c(4,0,1,1))
hist_vals <- table(tg_bins)
barplot(hist_vals, axes=FALSE, ylim=c(0+1.5, length(hist_vals)-1.5), xlim=c(0, max(hist_vals)), space=0, horiz=TRUE, border=NA, col="grey70", names="", xlab="nGenes", bty="l")
axis(side=1)
axis(side=2, at=c(-1,seq(from=1, to=length(g_axis_labs), by=5)-0.5,60), 
	labels=rep("", times=length(seq(from=1, to=length(g_axis_labs), by=5))+2))

# blank R corner
source("~/R-Scripts/Blank_plot.R")
blank_plot()

# expr through time
par(mar=c(0,1,1,4))
plot(pseudotime, cell_spliced_s, xlab="", ylab="", 
	xlim=c(min(pseudotime)+0.22, max(pseudotime)-0.22), 
	xaxt="n", yaxt="n", col=cell_color,
	ylim=c(0,150))
axis(side=4)
abline(v=p_breaks, lty=3)
points(pseudotime, cell_spliced_s,col=cell_color)

# spliced peaks Matching
par(mar=c(4,1,1,4))
image(x=p_breaks, y=g_breaks, z=t(log2(heat_data_e+1)), ylab="", xlab="pseudotime", col=down_col, yaxt="n")
axis(4, at=g_axis_labs[seq(from=1, to=length(g_axis_labs), by=5)], labels=g_axis_labs[seq(from=1, to=length(g_axis_labs), by=5)])

# trans gene Histogram
par(mar=c(4,1,1,0))
hist_vals <- table(sg_bins)
barplot(-hist_vals, axes=FALSE, ylim=c(0+1.5, length(hist_vals)-1.5), xlim=c(-max(hist_vals), 0), space=0, horiz=TRUE, border=NA, col="grey70", names="", xlab="nGenes", bty="l")
axis(side=1, at=seq(from=0, to=-100, by=-20), labels=seq(from=0, to=-100, by=-20)*-1)
axis(side=4, at=c(-1,seq(from=1, to=length(g_axis_labs), by=5)-0.5,60), 
	labels=rep("", times=length(seq(from=1, to=length(g_axis_labs), by=5))+2))

dev.off()

