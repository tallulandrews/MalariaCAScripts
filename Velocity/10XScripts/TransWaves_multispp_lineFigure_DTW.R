# Examining transcriptional waves using RNA velocity through IDC

# Colours
require("RColorBrewer")
require("Matrix")
require("gplots")
require("methods")

tab_color_remap <- read.table("Cluster_color_remap.csv", sep=",", stringsAsFactors=FALSE)
remap_colors <- function(old_colors) {
	for (i in 1:nrow(tab_color_remap)) {
        	old_colors[old_colors == tab_color_remap[i,1]] <- tab_color_remap[i,2]
	}
	return(old_colors);
}


## 10X
# Bergei
bergei_dat_10X <- readRDS("bergei_10X_velodata_wPseudo.rds")
bergei_velo_10X <- readRDS("bergei_10X_veloEst.rds")

tmp<- bergei_velo_10X$deltaE
tmp[tmp < 0] <- 0
bergei_UP <- tmp
bergei_cell_up <- Matrix::colSums(tmp)
bergei_up_s <- Matrix::colSums( t(apply(tmp, 1, function(a) {a<- a-min(a); a <- a/max(a);})) )
bergei_colors <- bergei_dat_10X$ann$cell_color
bergei_pseudo <- bergei_dat_10X$ann$pseudotime

# Knowlesi
knowlesi_dat_10X <- readRDS("new_pk_10X_velodata.rds")
knowlesi_velo_10X <- readRDS("new_pk_10X_veloEst.rds")

tmp<- knowlesi_velo_10X$deltaE
tmp[tmp < 0] <- 0
knowlesi_cell_up <- Matrix::colSums(tmp)
knowlesi_up_s <- Matrix::colSums( t(apply(tmp, 1, function(a) {a<- a-min(a); a <- a/max(a);})) )
knowlesi_colors <- knowlesi_dat_10X$ann$cell_color
knowlesi_pseudo <- knowlesi_dat_10X$ann$pseudotime
knowlesi_scmap_pseudo <- knowlesi_dat_10X$ann
knowl_cleanup <- knowlesi_colors != "black"
knowlesi_scmap_pseudo <- knowlesi_scmap_pseudo[knowl_cleanup,]
knowlesi_colors <- knowlesi_colors[knowl_cleanup]
knowlesi_cell_up <- knowlesi_cell_up[knowl_cleanup]
knowlesi_up_s <- knowlesi_up_s[knowl_cleanup]


# Falciparum

falciparum_dat_10X <- readRDS("new_pf_10X_velodata.rds")
falciparum_velo_10X <- readRDS("new_pf_10X_veloEst.rds")
#falciparum_dat_10X <- readRDS("noncomp_pf_10X_velodata.rds")
#falciparum_velo_10X <- readRDS("noncomp_pf_10X_veloEst.rds")

tmp<- falciparum_velo_10X$deltaE
tmp[tmp < 0] <- 0
falciparum_cell_up <- Matrix::colSums(tmp)
falciparum_up_s <- Matrix::colSums( t(apply(tmp, 1, function(a) {a<- a-min(a); a <- a/max(a);})) )
falciparum_colors <- falciparum_dat_10X$ann$cell_color
falciparum_pseudo <- falciparum_dat_10X$ann$pseudotime
falciparum_scmap_pseudo <- falciparum_dat_10X$ann
falcip_cleanup <- falciparum_colors != "black"
falciparum_scmap_pseudo <- falciparum_scmap_pseudo[falcip_cleanup,]
falciparum_colors <- falciparum_colors[falcip_cleanup]
falciparum_cell_up <- falciparum_cell_up[falcip_cleanup]
falciparum_up_s <- falciparum_up_s[falcip_cleanup]


# Bergei SCE
pb_mat <- bergei_dat_10X$emat+bergei_dat_10X$smat
pb_ngenes_c <- Matrix::colSums(pb_mat>0)
pb_numi_g_c <- Matrix::colSums(pb_mat)/pb_ngenes_c

# Knowlesi SCE
pk_mat <- knowlesi_dat_10X$emat+knowlesi_dat_10X$smat
pk_mat <- pk_mat[,knowl_cleanup]
pk_ngenes_c <- Matrix::colSums(pk_mat>0)
pk_numi_g_c <- Matrix::colSums(pk_mat)/pk_ngenes_c


# Falciparum SCE
pf_mat <- falciparum_dat_10X$emat+falciparum_dat_10X$smat
pf_mat <- pf_mat[,falcip_cleanup]
pf_ngenes_c <- Matrix::colSums(pf_mat>0)
pf_numi_g_c <- Matrix::colSums(pf_mat)/pf_ngenes_c

# Color Bars
# expression blocks

# Cell-specific
bergei_cell_breaks <- seq(from=min(bergei_pseudo), to=max(bergei_pseudo), length=151)
pseudo_cell_bins <- cut(bergei_pseudo, breaks=bergei_cell_breaks, include.lowest=TRUE)
pseudo_cell_bins_ord <- pseudo_cell_bins[order(bergei_pseudo)]
bar_pos <- cbind(bergei_cell_breaks[1:(length(bergei_cell_breaks)-1)],
		 rep(0, length(bergei_cell_breaks)-1),
		 bergei_cell_breaks[2:length(bergei_cell_breaks)],
		 rep(1, length(bergei_cell_breaks)-1)
		)
# painter
extern <- read.delim("/nfs/users/nfs_t/ta6/Collaborations/MCA/Pf_transcription_vs_time_FromSuppl1.csv", sep=",", header=T)
extern <- extern[!is.na(extern[,51]), ]
ortho <- read.table("/nfs/users/nfs_t/ta6/Collaborations/MCA/orthologueTable.csv", sep=",", header=T)
ortho <- ortho[ortho[,4] %in% extern$Gene.ID,]

bergei_UP <- bergei_UP[rownames(bergei_UP) %in% ortho$Gene,]
extern <- extern[extern$Gene.ID %in% ortho[,4],]
extern$pbanka <- ortho[match(extern$Gene.ID, ortho[,4]), 1]
extern <- extern[match(rownames(bergei_UP), extern$pbanka),]
bergei_UP_s <- t(apply(bergei_UP, 1, function(a) {a<- a-min(a); a <- a/max(a);}))
peak_time <- extern[,51]
#cell_peaks <- apply(bergei_UP_s, 2, function(a) {sum(a*peak_time)/sum(a)}) # diff num genes through time



n_gbins=45
g_breaks <- seq(from=min(extern[,51]), to=max(extern[,51]), length=n_gbins+1)
g_bins <- cut(extern[,51], breaks=n_gbins);
require("CellTypeProfiles")
data_s <- t(my_row_mean_aggregate(t(bergei_UP_s), g_bins));
bin_time <- (g_breaks[1:(length(g_breaks)-1)]+g_breaks[2:(length(g_breaks))] )/2
bergei_cell_peaks <- apply(data_s, 2, function(a) {sum(a*bin_time)/sum(a)})

cell_peak_ord <- bergei_cell_peaks[order(bergei_pseudo)]
cell_peak_bar_vals <- aggregate(cell_peak_ord, list(pseudo_cell_bins_ord), mean)
poor_match <- apply(data_s, 2, max)<0.1
poor_match_bin <- aggregate(poor_match[order(bergei_pseudo)], list(pseudo_cell_bins_ord), mean)

peak_cols <- colorRampPalette(c("mistyrose", "pink", "orchid1", "orchid3", "orchid4", "violetred1", "violetred3","violetred4"))(20)
peak_cols <- colorRampPalette(c("white", "#252525"))(8)
peak_bar_cols <- peak_cols[cut(cell_peak_bar_vals[,2], breaks=length(peak_cols))]
peak_bar_cols[poor_match_bin[,2]>=0.5] <- "cadetblue1"

# alternate time assignment - correlations
require("Hmisc")
out <- rcorr(as.matrix(log2(extern[,3:50]+1)), as.matrix(log2(bergei_UP+1)), type="pearson")
thing <- out$r[1:48,-c(1:48)]
assigned <- apply(thing,2, function(a){which(a==max(a))})
assigned_r <- apply(thing,2, function(a){a[which(a==max(a))]})

assigned_ord <- assigned[order(bergei_pseudo)]
assigned_bar_vals <- aggregate(assigned_ord, list(pseudo_cell_bins_ord), mean)
assigned_bar_cols <- peak_cols[cut(assigned_bar_vals[,2], breaks=length(peak_cols))]

## DTW ##
pseudotimes <- list(pb=bergei_pseudo, pk=knowlesi_scmap_pseudo$pseudotime, pf=falciparum_scmap_pseudo$pseudotime)

# Create matrices of only increasing expression profiles for each species
pb_incr <- bergei_velo_10X$deltaE
pb_incr[pb_incr<0] <- 0
pb_incr <- t(apply(pb_incr, 1, function(a) {a<- a-min(a); a <- a/max(a);}))
pk_incr <- knowlesi_velo_10X$deltaE[,knowl_cleanup]
pk_incr[pk_incr<0] <- 0
pk_incr <- t(apply(pk_incr, 1, function(a) {a<- a-min(a); a <- a/max(a);}))
pf_incr <- falciparum_velo_10X$deltaE[,falcip_cleanup]
pf_incr[pf_incr<0] <- 0
pf_incr <- t(apply(pf_incr, 1, function(a) {a<- a-min(a); a <- a/max(a);}))


# Smoothing
width <- 2*pi/100
bins <- seq(from=0, to=2*pi, length=500)

sm_fun <- function(b, incr, time, width){
	subset <- time > b-width/2 & time < b+width/2;
	if (sum(subset) > 1) {
		return(Matrix::rowMeans(incr[,subset]))
	} else if (sum(subset) == 1) {
		return(incr[,subset])
	} else {
		return(rep(NA, nrow(incr)));
	}
}
 
pb_incr_sm <- sapply(bins, sm_fun, incr=pb_incr, time=pseudotimes$pb, width=width)
rownames(pb_incr_sm) <- rownames(pb_incr)
pk_incr_sm <- sapply(bins, sm_fun, incr=pk_incr, time=pseudotimes$pk, width=width)
rownames(pk_incr_sm) <- rownames(pk_incr)
pf_incr_sm <- sapply(bins, sm_fun, incr=pf_incr, time=pseudotimes$pf, width=width)
rownames(pf_incr_sm) <- rownames(pf_incr)

time2bins <- function(t, bins, width) { return(which( (bins-width) < t & (bins+width) > t)) }
pb_to_bins <- sapply(pseudotimes$pb, time2bins, bins=bins, width=width)
pk_to_bins <- sapply(pseudotimes$pk, time2bins, bins=bins, width=width)
pf_to_bins <- sapply(pseudotimes$pf, time2bins, bins=bins, width=width)


# Orthology_match
#extern <- read.delim("/nfs/users/nfs_t/ta6/Collaborations/MCA/Pf_transcription_vs_time_FromSuppl1.csv", sep=",", header=T)
#extern <- extern[!is.na(extern[,51]), ]
#ortho <- read.table("/nfs/users/nfs_t/ta6/Collaborations/MCA/orthologueTable.csv", sep=",", header=T)
ortho <- read.table("/nfs/users/nfs_t/ta6/Collaborations/MCA/OrthologTable_20190413.csv", sep=",", header=T)
ortho$pbanka <- ortho[,2]
ortho$pfalci <- ortho[,4]
ortho$pknowl <- ortho[,3]
#ortho$pbanka <- ortho[,1]
#ortho$pfalci <- ortho[,4]
#ortho$pknowl <- ortho[,6]
#ortho$pknowl <- sub("PKH", "PKNH", ortho$pknowl)
#ortho$pknowl <- paste(ortho$pknowl, "0", sep="")

names(bins) <- paste("bin", 1:length(bins));
pk_sm_orth <- pk_incr_sm[match(ortho$pknowl, rownames(pk_incr_sm)),]
colnames(pk_sm_orth) <- names(bins);
pf_sm_orth <- pf_incr_sm[match(ortho$pfalci, rownames(pf_incr_sm)),]
colnames(pf_sm_orth) <- names(bins);
pb_sm_orth <- pb_incr_sm[match(ortho$pbanka, rownames(pb_incr_sm)),]
colnames(pb_sm_orth) <- names(bins);
pb_sm_to_pf <- pb_sm_orth[!is.na(rownames(pf_sm_orth)) & !is.na(rownames(pb_sm_orth)),]
pf_sm_to_pb <- pf_sm_orth[!is.na(rownames(pf_sm_orth)) & !is.na(rownames(pb_sm_orth)),]
pb_sm_to_pk <- pb_sm_orth[!is.na(rownames(pk_sm_orth)) & !is.na(rownames(pb_sm_orth)),]
pk_sm_to_pb <- pk_sm_orth[!is.na(rownames(pk_sm_orth)) & !is.na(rownames(pb_sm_orth)),]
pf_sm_to_pb <- pf_sm_to_pb[,!is.na(colSums(pf_sm_to_pb))]



# dtw by expression correlation-distance
require(dtw)
source("~/Collaborations/MCA/Velocity/10XScripts/plot.dtw.R")

# rows = query, columns = reference
d_bk <- function(i,j) {	return(1-cor(pk_sm_to_pb[,i], pb_sm_to_pk[,j], method="spearman")) }
d_bf <- function(i,j) {	return(1-cor(pf_sm_to_pb[,i], pb_sm_to_pf[,j], method="spearman")) }

# cheat
#d_bk <- function(i,j) { a <- colSums(pk_sm_to_pb); b <- colSums(pb_sm_to_pk); return(abs(a[i]/mean(a)-b[j]/mean(b)))}
#d_bf <- function(i,j) { a <- colSums(pf_sm_to_pb); b <- colSums(pb_sm_to_pf); return(abs(a[i]/mean(a)-b[j]/mean(b)))}


# b vs k
mat_bk <- matrix(0, nrow=ncol(pk_sm_to_pb), ncol=ncol(pb_sm_to_pk));
for (i in 1:nrow(mat_bk)) {
	for(j in 1:ncol(mat_bk)) {
		mat_bk[i,j] <- d_bk(i,j)
	}
}

#pdf("Pk_to_Pb_DTW_heatmap.pdf", width=6, height=6)
heatmap.2(mat_bk, Rowv=FALSE, Colv=FALSE, trace="none", dendrogram="none", scale="column", col=colorRampPalette(c("grey90", "black"))(30), key.title="", key.xlab="Distance", xlab="P. bergei (t)", ylab="P. knowlesi (t)")
#dev.off()

test_bk <- dtw(mat_bk, keep.internals=TRUE) 
test_bk$query <- colSums(pk_sm_to_pb)
test_bk$reference <- colSums(pb_sm_to_pk)

#pdf("Pk_to_Pb_DTW.pdf", width=4, height=4)
my_dtwPlotTwoWay(test_bk, xlab="Time", ylab="Mean increase in expression", main="Pk vs Pb")
#dev.off()

png("Pk_to_Pb_DTW_densityContour.png", width=6, height=6, units="in", res=300)
my_dtwPlotDensity(test_bk, normalize=TRUE, xlab="Pk Time Index", ylab="Pb Time Index")
#dev.off()


# index 1 = query index
# index 2 = reference index
# b vs f
mat_bf <- matrix(0, nrow=ncol(pf_sm_to_pb), ncol=ncol(pb_sm_to_pf));
for (i in 1:nrow(mat_bf)) {
	for(j in 1:ncol(mat_bf)) {
		mat_bf[i,j] <- d_bf(i,j)
	}
}

#pdf("Pf_to_Pb_DTW_heatmap.pdf", width=6, height=6)
heatmap.2(mat_bf, Rowv=FALSE, Colv=FALSE, trace="none", dendrogram="none", scale="column", col=colorRampPalette(c("grey90", "black"))(30), key.title="", key.xlab="Distance", xlab="P. bergei (t)", ylab="P. falciparum (t)")
#dev.off()

test_bf <- dtw(mat_bf, keep.internals=TRUE) 
test_bf$query <- colSums(pf_sm_to_pb)
test_bf$reference <- colSums(pb_sm_to_pf)
deltaQs <- vector(length=length(unique(test_bf$index1)));
names(deltaQs) <- colnames(pf_sm_to_pb)
for (i in unique(test_bf$index1)) {
	Refs <- unique(test_bf$index2[test_bf$index1==i])
	deltaQ <- mean(bins[Refs]) - bins[i]
	deltaQs[i] <- deltaQ
}

pf_shifts <- lapply(pf_to_bins, function(b, bins){mean(deltaQs[names(bins)[b]], na.rm=TRUE)}, bins=bins)

pf_pseudotime_shifted <- pseudotimes$pf+unlist(pf_shifts)
# Just in case:
pf_pseudotime_shifted[pf_pseudotime_shifted < 0] <- 0;
pf_pseudotime_shifted[pf_pseudotime_shifted > 2*pi] <- 2*pi;

#pdf("Pf_to_Pb_DTW.pdf", width=4, height=4)
my_dtwPlotTwoWay(test_bf,scalex=3.5, xlab="Time", ylab="Mean increase in expression", main="Pf vs Pb")
#dev.off()

png("Pf_to_Pb_DTW_densityContour.png", width=6, height=6, units="in", res=300)
my_dtwPlotDensity(test_bf, normalize=TRUE, xlab="Pf Time Index", ylab="Pb Time Index")
#dev.off()


#plot(colSums(pf_sm_to_pb)[test_bf$index2])
#plot(colSums(pb_sm_to_pf)[test_bf$index1])

#plot(pf_pseudotime_shifted, colSums(pf_incr))
#plot(pseudotimes$pb, colSums(pb_incr))


#pseudotimes$pf <- pf_pseudotime_shifted;




### Ngene/umi ###
# ngene
ngene_cols <- colorRampPalette(c("white", "#252525"))(8)
ngene_vals <- pb_ngenes_c
ngene_ord <- ngene_vals[order(pseudotimes$pb)]
ngene_bar_vals <- aggregate(ngene_ord, list(pseudo_cell_bins_ord), mean)
ngene_bar_cols <- ngene_cols[cut(ngene_bar_vals[,2], breaks=length(ngene_cols))]


# numi/gene
numi_cols <- colorRampPalette(c("white", "#252525"))(8)
numi_vals <- pb_numi_g_c
numi_ord <- numi_vals[order(pseudotimes$pb)]
numi_bar_vals <- aggregate(numi_ord, list(pseudo_cell_bins_ord), mean)
numi_bar_cols <- numi_cols[cut(numi_bar_vals[,2], breaks=length(numi_cols))]

# Colour scale bars
source("~/R-Scripts/Colour_bar.R")
#pdf("DTW_TransWaves_multispp_lineFigure_colorscales.pdf", width=2, height=6)
par(mfrow=c(1,3));
color.bar(peak_cols, ticks.at=seq(0,1,len=length(peak_cols)+1), 
		ticks.lab=round(seq(min(cell_peak_bar_vals[,2]), max(cell_peak_bar_vals[,2]), len=length(peak_cols)+1)), 
		title="", add=TRUE)
title(ylab="Painter(h)", line=2.5, cex.lab=1.5)

color.bar(ngene_cols, ticks.at=seq(0,1,len=length(ngene_cols)+1), 
		ticks.lab=round(seq(min(ngene_bar_vals[,2]), max(ngene_bar_vals[,2]), len=length(ngene_cols)+1)), 
		title="", add=TRUE)
title(ylab="nGene", line=2.5, cex.lab=1.5)

color.bar(numi_cols, ticks.at=seq(0,1,len=length(numi_cols)+1), 
		ticks.lab=round(seq(min(numi_bar_vals[,2]), max(numi_bar_vals[,2]), len=length(numi_cols)+1), digits=1), 
		title="", add=TRUE)
title(ylab="nUMI/g", line=2.5, cex.lab=1.5)

#dev.off()

plot_bar <- function(bar.colors, bar.pos, xlim=c(0,1)) { 
	plot(1,1, col="white", xlim=c(0,6.28), ylim=c(0,1), xaxt="n", yaxt="n", main="", xlab="", ylab="", bty="n")
	min_x <- min(bar.pos[,1], bar.pos[,3])
	max_x <- max(bar.pos[,1], bar.pos[,3])
	min_y <- min(bar.pos[,2], bar.pos[,4])
	max_y <- max(bar.pos[,2], bar.pos[,4])
	range <- xlim[2]-xlim[1]
	bar.pos[,1] <- (bar.pos[,1]- min_x)/(max_x-min_x)*range + xlim[1]
	bar.pos[,3] <- (bar.pos[,3]- min_x)/(max_x-min_x)*range + xlim[1]
	bar.pos[,2] <- (bar.pos[,2]- min_y)/(max_y-min_y)
	bar.pos[,4] <- (bar.pos[,4]- min_y)/(max_y-min_y)

	for (i in 1:length(bar.colors)) {
		rect(bar.pos[i,1], bar.pos[i,2], bar.pos[i,3], bar.pos[i,4], col=bar.colors[i], border=NA)
	}
}
		

# bin cells by pseudotime
n_pbins=20
p_breaks <- seq(from=min(pseudotimes$pb), to=max(pseudotimes$pb), length=n_pbins+1)
p_bins <- cut(pseudotimes$pb, breaks=n_pbins)

nlines=1
nspp=3
my_xlim <- c(min(p_breaks), max(p_breaks))

### Reviewer Response ###
# Analytic comparison of waves!

# Automatically split curve into segments
pb_smoothed <- smooth.spline(pseudotimes$pb, bergei_up_s)
bit <- pseudotimes$pk>-1
pk_smoothed <- smooth.spline(pseudotimes$pk[bit], knowlesi_up_s[bit], spar=0.9)
bit <- pseudotimes$pf > 0.5 & pseudotimes$pf < 6.5
pf_smoothed <- smooth.spline(pseudotimes$pf[bit], falciparum_up_s[bit], spar=0.9)

delta_y <- pb_smoothed$y[2:length(pb_smoothed$y)] - pb_smoothed$y[1:(length(pb_smoothed$y)-1)]
splits <- c(1);
for (i in 1:(length(delta_y)-1)) {
	if (sign(delta_y[i]) != sign(delta_y[i+1])){ splits <- c(splits, i)}
}
splits <- c(splits, length(delta_y));
splits <- pb_smoothed$x[splits]

plot_piecewise <- function(raw_x, raw_y, splits) {
	for (j in 1:(length(splits)-1)) {
		dy <- raw_y[raw_x > splits[j] & raw_x < splits[j+1]]
		dx <- raw_x[raw_x > splits[j] & raw_x < splits[j+1]]
		if (length(dy) < 3 | var(dx) == 0) {next;}
		reg <- lm(dy~dx)
		out <- summary(reg)
		if(out$coefficients[2,4] < 0.05/length(splits)){
			sty <- min(dx)*out$coefficients[2,1]+out$coefficients[1,1]
			endy <- max(dx)*out$coefficients[2,1]+out$coefficients[1,1]
			lines(c(min(dx), max(dx)), c(sty, endy), col="red", lwd=2, lty=2)
		}
	}
}

# Add Sig features to plot
nlines=1
pdf("RR_Three_Spp_Wave_Lines_w_peaktrough.pdf", width=7, height=9)
lmat <- cbind(c(1:nlines, nlines+1:nspp), c(1:nlines, nlines+1:nspp))
layout(lmat, width=c(1,1), heights=c(rep(0.1, times=nlines), rep(1, times=nspp)))

#lines
par(mar=c(0,4,0,1));
plot_bar(assigned_bar_cols, bar_pos, xlim=my_xlim)
mtext("Painter", side=2, line=-1.5, las=2)

#spp
# BERGEI
par(mar=c(0,4,0.2,1))
plot(pseudotimes$pb, bergei_up_s, col=remap_colors(bergei_colors), ylab="Bergei", xlab="", xaxt="n", xlim=my_xlim)
#abline(v=p_breaks, lty=3, col="grey40")
abline(v=splits, lty=3, col="grey40")
points(pseudotimes$pb, bergei_up_s, col=remap_colors(bergei_colors))
thing <- smooth.spline(pseudotimes$pb, bergei_up_s)
lines(thing$x, thing$y, lwd=2)

# sig
#plot_piecewise(pseudotimes$pb, bergei_up_s, splits)

# KNOWLESI
par(mar=c(0,4,0,1))
plot(pseudotimes$pk, knowlesi_up_s, col=remap_colors(knowlesi_colors), ylab="Knowlesi", xlab="", xaxt="n", xlim=my_xlim)
#abline(v=p_breaks, lty=3, col="grey40")
abline(v=splits, lty=3, col="grey40")
points(pseudotimes$pk, knowlesi_up_s, col=remap_colors(knowlesi_colors))
bit1 <- pseudotimes$pk>-1
thing1 <- smooth.spline(pseudotimes$pk[bit1], knowlesi_up_s[bit1], spar=0.9)
lines(thing1$x, thing1$y, lwd=2)
# sig
#plot_piecewise(pseudotimes$pk, knowlesi_up_s, splits)

#FALCIPARUM
par(mar=c(4,4,0,1))
plot(pseudotimes$pf, falciparum_up_s, col=remap_colors(falciparum_colors), xlab="Pseudotime", ylab="Falciparum", xlim=my_xlim)
#abline(v=p_breaks, lty=3, col="grey40")
abline(v=splits, lty=3, col="grey40")
points(pseudotimes$pf, falciparum_up_s, col=remap_colors(falciparum_colors))
bit1 <- pseudotimes$pf > 0.5 & pseudotimes$pf < 6.5
thing1 <- smooth.spline(pseudotimes$pf[bit1], falciparum_up_s[bit1], spar=0.9)
lines(thing1$x, thing1$y, lwd=2)
# sig
#plot_piecewise(pseudotimes$pf, falciparum_up_s, splits)

dev.off()


# DE genes for pieces are the same across species?
get_de_mat <- function(g_ortho, this_incr, this_time) {
	MAT <- matrix(0, nrow=length(g_ortho), ncol=length(splits)-1);
	MAT_p <- matrix(0, nrow=length(g_ortho), ncol=length(splits)-1);
	MAT_v <- matrix(0, nrow=length(g_ortho), ncol=length(splits)-1);
	rownames(MAT) <- g_ortho;
	rownames(MAT_p) <- g_ortho;
	rownames(MAT_v) <- g_ortho;
	for (i in 2:length(splits)) {
		test_gene <- function(g,i, incr, time) {
			yes <-incr[rownames(incr)==g, time < splits[i] & time > splits[i-1]]
			xes <- time[time < splits[i] & time > splits[i-1]]
			out <- summary(lm(yes~xes))$coefficients[2,c(1,4)]
			junk <- yes[1:2]; junk[1]=0; junk[2]=1;
			if (is.na(out[1]) | is.na(out[2])) {return(junk)}
			return(out)
		}
		out<-t(sapply(g_ortho, test_gene, i = i, incr=this_incr, time=this_time))
		out <- as.data.frame(out)
		out$q.value <- p.adjust(out[,2], method="fdr")
		sig <- out[out$q.value < 0.05,]
		MAT[rownames(MAT) %in% rownames(sig),i-1] <- sig$Estimate
		MAT_v[,i-1] <- out$Estimate
		MAT_p[,i-1] <- out$q.value
	}
	return(list(sig=MAT, est=MAT_v, qval=MAT_p));
}
g_ortho <- union(rownames(pb_sm_to_pk), rownames(pb_sm_to_pf))
MAT_pb <- get_de_mat(rownames(pb_incr), pb_incr, pseudotimes$pb)
saveRDS(MAT_pb, file="pb_Velo_MAT_de.rds")
g_ortho <- rownames(pk_sm_to_pb)
MAT_pk <- get_de_mat(rownames(pk_incr), pk_incr, pseudotimes$pk)
saveRDS(MAT_pk, file="pk_Velo_MAT_de.rds")
g_ortho <- rownames(pf_sm_to_pb)
MAT_pf <- get_de_mat(rownames(pf_incr), pf_incr, pseudotimes$pf)
saveRDS(MAT_pf, file="pf_Velo_MAT_de.rds")
print("PiecewiseDE Done")

MAT_pb <- MAT_pb$sig
MAT_pk <- MAT_pk$sig
MAT_pf <- MAT_pf$sig

ortho <- read.table("/nfs/users/nfs_t/ta6/Collaborations/MCA/OrthologTable_20190413.csv", sep=",", header=T)
ortho$pbanka <- ortho[,2]
ortho$pfalci <- ortho[,4]
ortho$pknowl <- ortho[,3]

rownames(MAT_pk) <- as.character(ortho$pbanka[match(rownames(MAT_pk), ortho$pknowl)])
rownames(MAT_pf) <- as.character(ortho$pbanka[match(rownames(MAT_pf), ortho$pfalci)])

# Order
tmp <- MAT_pb;
for (i in seq(1,10,2)) {
	tmp[tmp[,i] < 0,i] <- 0
}
for (i in seq(2,10,2)) {
	tmp[tmp[,i] > 0,i] <- 0
}
tmp <- tmp[,1:10]
assign <- apply(tmp, 1, function(x){a <- which(abs(x)==max(abs(x))); if (length(a) > 1) {return(11)}else{return(a)}})
MAT_pb <- MAT_pb[order(assign),]
heatcol <- colorRampPalette(c("blue", "white", "red"))(30)
heatmap.2(MAT_pb, Colv=FALSE, Rowv=FALSE, trace="none", col=heatcol)

MAT_pf <- MAT_pf[match(rownames(MAT_pb), rownames(MAT_pf)),]
heatmap.2(MAT_pf, Colv=FALSE, Rowv=FALSE, trace="none", col=heatcol)
MAT_pk <- MAT_pk[match(rownames(MAT_pb), rownames(MAT_pk)),]
heatmap.2(MAT_pk, Colv=FALSE, Rowv=FALSE, trace="none", col=heatcol)


## PB vs PK
col_scale <- c(-100, seq(from=-.5, to=.5, length=29), 100)
pb <- MAT_pb[rownames(MAT_pb) %in% rownames(MAT_pk),]
pk <- MAT_pk[match(rownames(pb), rownames(MAT_pk)),]
png("piecewise_pb_pk_de_part1.png", width=5, height=5, units="in", res=300)
heatmap.2(pb, Colv=FALSE, Rowv=FALSE, trace="none", col=heatcol, breaks=col_scale)
#dev.off()
#png("piecewise_pb_pk_de_part2.png", width=5, height=5, units="in", res=300)
heatmap.2(pk, Colv=FALSE, Rowv=FALSE, trace="none", col=heatcol, breaks=col_scale)
#dev.off()

## PB vs PF
pb <- MAT_pb[rownames(MAT_pb) %in% rownames(MAT_pf),]
pf <- MAT_pf[match(rownames(pb), rownames(MAT_pf)),]
#png("piecewise_pb_pf_de_part1.png", width=5, height=5, units="in", res=300)
heatmap.2(pb, Colv=FALSE, Rowv=FALSE, trace="none", col=heatcol, breaks=col_scale)
#dev.off()
#png("piecewise_pb_pf_de_part2.png", width=5, height=5, units="in", res=300)
heatmap.2(pf, Colv=FALSE, Rowv=FALSE, trace="none", col=heatcol, breaks=col_scale)
#dev.off()

### Triple Figure ###
consistent_genes <- intersect(as.character(rownames(pb)), intersect(as.character(rownames(pf)), as.character(rownames(pk))))

pb_trip <- pb[rownames(pb) %in% consistent_genes,]
pk_trip <- pk[rownames(pk) %in% consistent_genes,]
pf_trip <- pf[rownames(pf) %in% consistent_genes,]

#pdf("Triple_Align_cycle_heatmap.pdf", width=11, height=8)
triple_xlim=c(0+0.23, 2*pi-0.25)
x_loc <- splits[1:(length(splits)-1)]+diff(splits)
x_lab <- 1:(length(splits)-1)
layout(rbind(c(1,2,3), c(4,5,6)), widths=c(1,1,1), heights=c(1,4))

par(mar=c(0, 1, 2, 1))
# pb line
plot(pseudotimes$pb, bergei_up_s, col=remap_colors(bergei_colors), main="Bergei", xlab="", ylab="", xaxt="n", yaxt="n", xlim=triple_xlim)
abline(v=splits, lty=3, col="grey40")
points(pseudotimes$pb, bergei_up_s, col=remap_colors(bergei_colors))
thing <- smooth.spline(pseudotimes$pb, bergei_up_s)
lines(thing$x, thing$y, lwd=2)
plot_piecewise(pseudotimes$pb, bergei_up_s, splits)

# pk line
plot(pseudotimes$pk, knowlesi_up_s, col=remap_colors(knowlesi_colors), main="Knowlesi", xlab="", ylab="", xaxt="n", yaxt="n", xlim=triple_xlim)
#abline(v=p_breaks, lty=3, col="grey40")
abline(v=splits, lty=3, col="grey40")
points(pseudotimes$pk, knowlesi_up_s, col=remap_colors(knowlesi_colors))
bit1 <- pseudotimes$pk>-1
thing1 <- smooth.spline(pseudotimes$pk[bit1], knowlesi_up_s[bit1], spar=0.9)
lines(thing1$x, thing1$y, lwd=2)
plot_piecewise(pseudotimes$pk, knowlesi_up_s, splits)

# pf line
plot(pseudotimes$pf, falciparum_up_s, col=remap_colors(falciparum_colors), xlab="Pseudotime", main="Falciparum", ylab="", xlim=triple_xlim, xaxt="n", yaxt="n")
#abline(v=p_breaks, lty=3, col="grey40")
abline(v=splits, lty=3, col="grey40")
points(pseudotimes$pf, falciparum_up_s, col=remap_colors(falciparum_colors))
bit1 <- pseudotimes$pf > 0.5 & pseudotimes$pf < 6.5
thing1 <- smooth.spline(pseudotimes$pf[bit1], falciparum_up_s[bit1], spar=0.9)
lines(thing1$x, thing1$y, lwd=2)
plot_piecewise(pseudotimes$pf, falciparum_up_s, splits)

par(mar=c(4, 1, 0, 1))
# pb heatmap
image(x=splits,y = 1:nrow(pb_trip), t(pb_trip)[,rev(1:nrow(pb_trip))], col=heatcol, breaks=col_scale, xlab="Pseudotime", ylab="Genes", yaxt="n")
abline(v=splits, lty=3, col="grey40")

# pk heatmap
image(x=splits,y = 1:nrow(pk_trip), t(pk_trip)[,rev(1:nrow(pk_trip))], col=heatcol, breaks=col_scale, xlab="Pseudotime", ylab="Genes", yaxt="n")
abline(v=splits, lty=3, col="grey40")

# pf heatmap
image(x=splits,y = 1:nrow(pf_trip), t(pf_trip)[,rev(1:nrow(pf_trip))], col=heatcol, breaks=col_scale, xlab="Pseudotime", ylab="Genes", yaxt="n")
abline(v=splits, lty=3, col="grey40")

#dev.off()
# stats
for (i in 1:ncol(pb_trip)) {
s <- sum(sign(pb_trip[,i])*sign(pf_trip[,9]) == -1); o <- s/sum(pb_trip[,i] != 0)
print(c(i, o, s))}

source("~/R-Scripts/Colour_bar.R")
#pdf("triple_heatmap_colour_scale.pdf", width=1.5, height=3)
color.bar(heatcol, min=-0.6, max=0.6, ticks.at=seq(-0.6, 0.6, 0.15), add=TRUE)
title(ylab="piecewise slope")
#dev.off()



### Tables ###
ortho <- read.table("/nfs/users/nfs_t/ta6/Collaborations/MCA/OrthologTable_20190413.csv", sep=",", header=T)
ortho$pbanka <- ortho[,2]
ortho$pfalci <- ortho[,4]
ortho$pknowl <- ortho[,3]

time_mid <- signif(splits[1:(length(splits)-1)]+diff(splits)/2, 2)

MAT_pb<-readRDS("pb_Velo_MAT_de.rds")
MAT_pk<-readRDS("pk_Velo_MAT_de.rds")
MAT_pf<-readRDS("pf_Velo_MAT_de.rds")



prettify <- function (MAT) {
	out <- matrix(0, nrow=nrow(MAT$sig), ncol=2*ncol(MAT$sig))
	rownames(out) <- rownames(MAT$est);
	colnames(MAT$est) <- paste(as.character(time_mid), "slope", sep="_");
	colnames(MAT$qval) <- paste(as.character(time_mid), "qval", sep="_");

	out[,seq(1,ncol(out), 2)] <- MAT$est;
	out[,seq(2,ncol(out), 2)] <- MAT$qval;
	colnames(out) <- as.vector(rbind(colnames(MAT$est), colnames(MAT$qval)))
#	tmp <- paste(signif(MAT$est,2), signif(MAT$qval,2));
#	tmp <- sub(" ", " (", tmp);
#	tmp <- paste(tmp, ")", sep="");
#	new_mat <- matrix(tmp, ncol=ncol(MAT$est), byrow=FALSE);
#	rownames(new_mat) <- rownames(MAT$est)
	return(out);
}
pb_pretty <- prettify(MAT_pb)
#colnames(pb_pretty) <- as.character(time_mid)
pb_pretty <- data.frame(pb_pretty);
pb_pretty$pf <- ortho[match(rownames(pb_pretty), ortho$pbanka), "pfalci"]
pb_pretty$pk <- ortho[match(rownames(pb_pretty), ortho$pbanka), "pknowl"]
#pb_pretty$in_heatmap <- rownames(pb_pretty) %in% consistent_genes
write.table(pb_pretty, file="RR_full_pb_table.csv", sep=",")

pk_pretty <- prettify(MAT_pk)
#colnames(pk_pretty) <- as.character(time_mid)
pk_pretty <- data.frame(pk_pretty);
pk_pretty$pf <- ortho[match(rownames(pk_pretty), ortho$pknowl), "pfalci"]
pk_pretty$pb <- ortho[match(rownames(pk_pretty), ortho$pknowl), "pbanka"]
write.table(pk_pretty, file="RR_full_pk_table.csv", sep=",")

pf_pretty <- prettify(MAT_pf)
#colnames(pf_pretty) <- as.character(time_mid)
pf_pretty <- data.frame(pf_pretty);
pf_pretty$pk <- ortho[match(rownames(pf_pretty), ortho$pfalci), "pknowl"]
pf_pretty$pb <- ortho[match(rownames(pf_pretty), ortho$pfalci), "pbanka"]
write.table(pf_pretty, file="RR_full_pf_table.csv", sep=",")

pb_pretty <- pb_pretty[match(rownames(pb_trip), rownames(pb_pretty)),]
pk_pretty <- pk_pretty[match(pb_pretty$pk, rownames(pk_pretty)),]
pf_pretty <- pf_pretty[match(pb_pretty$pf, rownames(pf_pretty)),]
mega_table <- vector();
mega_colnames <- vector();
for (i in 1:length(time_mid)) {
	fst <- (i-1)*2+1
	snd <- (i-1)*2+2
	mega_table <- cbind(mega_table, pb_pretty[,fst], pb_pretty[,snd], 
					pk_pretty[,fst], pk_pretty[,snd], 
					pf_pretty[,fst], pf_pretty[,snd])
	mega_colnames <- c(mega_colnames, 
		paste("pb", colnames(pb_pretty)[c(fst,snd)], sep="_"), 
		paste("pk", colnames(pk_pretty)[c(fst,snd)], sep="_"),
		paste("pf", colnames(pf_pretty)[c(fst,snd)], sep="_"))
	
	agree_bk <- sign(pb_trip[,i])*sign(pk_trip[,i]) == 1
	agree_bf <- sign(pb_trip[,i])*sign(pf_trip[,i]) == 1
	mega_table <- cbind(mega_table, agree_bk, agree_bf);
	mega_colnames <- c(mega_colnames, 
		paste(time_mid[i],"Pb_Pk_agree", sep="_"), 
		paste(time_mid[i], "Pb_Pf_agree", sep="_"))
}
colnames(mega_table) <- mega_colnames
mega_table <- data.frame(mega_table);
mega_table$pk <- rownames(pk_pretty)
mega_table$pf <- rownames(pf_pretty)
write.table(mega_table, file="RR_mega_table_heatmap.csv", sep=",")


## Olaps?

jaccard <- function(a, b){length(intersect(a,b))/ length(union(a,b))}

mat_up_bk <- matrix(0, nrow=ncol(MAT_pb), ncol=ncol(MAT_pk));
for (i in 1:nrow(mat_up_bk)) {
	for(j in 1:ncol(mat_up_bk)) {
		mat_up_bk[i,j] <- jaccard(rownames(MAT_pb)[MAT_pb[,i] > 0], rownames(MAT_pk)[MAT_pk[,j] > 0])
	}
}

mat_down_bk <- matrix(0, nrow=ncol(MAT_pb), ncol=ncol(MAT_pk));
for (i in 1:nrow(mat_down_bk)) {
	for(j in 1:ncol(mat_down_bk)) {
		mat_down_bk[i,j] <- jaccard(rownames(MAT_pb)[MAT_pb[,i] < 0], rownames(MAT_pk)[MAT_pk[,j] < 0])
	}
}

mat_up_bf <- matrix(0, nrow=ncol(MAT_pb), ncol=ncol(MAT_pf));
for (i in 1:nrow(mat_up_bf)) {
	for(j in 1:ncol(mat_up_bf)) {
		mat_up_bf[i,j] <- jaccard(rownames(MAT_pb)[MAT_pb[,i] > 0], rownames(MAT_pf)[MAT_pf[,j] > 0])
	}
}

mat_down_bf <- matrix(0, nrow=ncol(MAT_pb), ncol=ncol(MAT_pf));
for (i in 1:nrow(mat_down_bf)) {
	for(j in 1:ncol(mat_down_bf)) {
		mat_down_bf[i,j] <- jaccard(rownames(MAT_pb)[MAT_pb[,i] < 0], rownames(MAT_pf)[MAT_pf[,j] < 0])
	}
}

mat_both_bk <- mat_up_bk; mat_both_bk[,seq(2,10,2)] <- -1*mat_down_bk[,seq(2,10,2)]

key_genes <- c(
	rownames(MAT_pb)[MAT_pb[,1] > 0],
	rownames(MAT_pb)[MAT_pb[,2] < 0],
	rownames(MAT_pb)[MAT_pb[,3] > 0],
	rownames(MAT_pb)[MAT_pb[,4] < 0],
	rownames(MAT_pb)[MAT_pb[,5] > 0],
	rownames(MAT_pb)[MAT_pb[,6] < 0],
	rownames(MAT_pb)[MAT_pb[,7] > 0],
	rownames(MAT_pb)[MAT_pb[,8] < 0],
	rownames(MAT_pb)[MAT_pb[,9] > 0],
	rownames(MAT_pb)[MAT_pb[,10] < 0]
	)

pb_to_pk_heatmap <- pb_incr[match(rownames(pb_sm_to_pk), rownames(pb_incr)),order(pseudotimes$pb)]
pb_to_pf_heatmap <- pb_incr[match(rownames(pb_sm_to_pf), rownames(pb_incr)),order(pseudotimes$pb)]

pk_to_pb_heatmap <- pk_incr[match(rownames(pk_sm_to_pb), rownames(pk_incr)),order(pseudotimes$pk)]
pf_to_pb_heatmap <- pf_incr[match(rownames(pf_sm_to_pb), rownames(pf_incr)),order(pseudotimes$pf)]

score <-( MAT_bp %*% 1:13 ) 



## Pb vs Pk ##
#par(mfrow=c(1,2))
require("gplots")
a <- heatmap.2(pb_to_pk_heatmap[rownames(pb_to_pk_heatmap) %in% rownames(MAT_bp[MAT_bp[,1]!=0,]),], Colv=FALSE, Rowv=FALSE, trace="none", dendrogram="none")
a <- heatmap.2(pk_to_pb_heatmap[rownames(pb_to_pk_heatmap) %in% rownames(MAT_bp[MAT_bp[,1]!=0,]),], Colv=FALSE, Rowv=FALSE, trace="none", dendrogram="none")


