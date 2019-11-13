# Examining transcriptional waves using RNA velocity through IDC

# Colours
require("RColorBrewer")
require("Matrix")
require("gplots")

tab_color_remap <- read.table("Cluster_color_remap.csv", sep=",", stringsAsFactors=FALSE)
remap_colors <- function(old_colors) {
	for (i in 1:nrow(tab_color_remap)) {
        	old_colors[old_colors == tab_color_remap[i,1]] <- tab_color_remap[i,2]
	}
	return(old_colors);
}


## 10X
# gene_thresholds: detected > 1% of cells, max-min > 10;
count_qc <- function(count_mat) {
	lo <- apply(count_mat, 1, min);
	hi <- apply(count_mat, 1, max);
	perc <- Matrix::rowMeans(count_mat > 0)
	keep <- hi-lo > 3 & perc > 0.01
	return(count_mat[keep,]);
}
# Bergei
bergei_dat_10X <- readRDS("bergei_10X_velodata_wPseudo.rds")
bergei_velo_10X <- readRDS("bergei_10X_veloEst.rds")
# count
bergei_counts <- readRDS("pb_count_mat.rds")
colnames(bergei_counts) <- sub("-1$", "", colnames(bergei_counts))
bergei_counts <- bergei_counts[, match(colnames(bergei_dat_10X$emat), colnames(bergei_counts))]
bergei_counts <- count_qc(bergei_counts)


tmp<- bergei_velo_10X$deltaE
tmp[tmp < 0] <- 0
bergei_cell_up <- Matrix::colSums(tmp)
bergei_up_s <- Matrix::colSums( t(apply(tmp, 1, function(a) {a<- a-min(a); a <- a/max(a);})) )
bergei_colors <- bergei_dat_10X$ann$cell_color
bergei_pseudo <- bergei_dat_10X$ann$pseudotime
# velo
bergei_incr <- tmp
bergei_incr <- t(apply(bergei_incr, 1, function(a) {a<- a-min(a); a <- a/max(a);}))

# Knowlesi
knowlesi_dat_10X <- readRDS("new_pk_10X_velodata.rds")
knowlesi_velo_10X <- readRDS("new_pk_10X_veloEst.rds")
knowl_cleanup <- knowlesi_dat_10X$ann$cell_color != "black"
knowlesi_ann <- knowlesi_dat_10X$ann[knowl_cleanup,]
knowlesi_colors <- knowlesi_ann$cell_color

knowlesi_counts <- readRDS("new_pk_count_mat.rds")
colnames(knowlesi_counts) <- sub("-1$", "", colnames(knowlesi_counts))
knowlesi_counts <- knowlesi_counts[, match(rownames(knowlesi_ann), colnames(knowlesi_counts))]
knowlesi_counts <- count_qc(knowlesi_counts)

tmp<- knowlesi_velo_10X$deltaE[,knowl_cleanup]
tmp[tmp < 0] <- 0
knowlesi_cell_up <- Matrix::colSums(tmp)
knowlesi_up_s <- Matrix::colSums( t(apply(tmp, 1, function(a) {a<- a-min(a); a <- a/max(a);})) )
# velo
knowlesi_incr <- tmp
knowlesi_incr <- t(apply(knowlesi_incr, 1, function(a) {a<- a-min(a); a <- a/max(a);}))


# Falciparum
falciparum_dat_10X <- readRDS("new_pf_10X_velodata.rds")
falciparum_velo_10X <- readRDS("new_pf_10X_veloEst.rds")
#falciparum_dat_10X <- readRDS("noncomp_pf_10X_velodata.rds")
#falciparum_velo_10X <- readRDS("noncomp_pf_10X_veloEst.rds")
falcip_cleanup <- falciparum_dat_10X$ann$cell_color != "black"
falciparum_ann <- falciparum_dat_10X$ann[falcip_cleanup,]
falciparum_colors <- falciparum_ann$cell_color

falciparum_counts <- readRDS("new_pf_count_mat.rds")
colnames(falciparum_counts) <- sub("-1$", "", colnames(falciparum_counts))
falciparum_counts <- falciparum_counts[, match(rownames(falciparum_ann), colnames(falciparum_counts))]
falciparum_counts <- count_qc(falciparum_counts)

tmp<- falciparum_velo_10X$deltaE[,falcip_cleanup]
tmp[tmp < 0] <- 0
falciparum_cell_up <- Matrix::colSums(tmp)
falciparum_up_s <- Matrix::colSums( t(apply(tmp, 1, function(a) {a<- a-min(a); a <- a/max(a);})) )
# velo
falciparum_incr <- tmp
falciparum_incr <- t(apply(falciparum_incr, 1, function(a) {a<- a-min(a); a <- a/max(a);}))

dat_version="Expr"
if (dat_version=="Velo") {
	bergei_counts <- bergei_incr
	knowlesi_counts <- knowlesi_incr
	falciparum_counts <- falciparum_incr
}

## DTW ##
pseudotimes <- list(pb=bergei_pseudo, pk=knowlesi_ann$pseudotime, pf=falciparum_ann$pseudotime)
trans_rates <- list(pb=bergei_up_s, pk=knowlesi_up_s, pf=falciparum_up_s);

# Create matrices of only increasing expression profiles for each species



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
 
pb_incr_sm <- sapply(bins, sm_fun, incr=bergei_counts, time=pseudotimes$pb, width=width)
rownames(pb_incr_sm) <- rownames(bergei_counts)
pk_incr_sm <- sapply(bins, sm_fun, incr=knowlesi_counts, time=pseudotimes$pk, width=width)
rownames(pk_incr_sm) <- rownames(knowlesi_counts)
pf_incr_sm <- sapply(bins, sm_fun, incr=falciparum_counts, time=pseudotimes$pf, width=width)
rownames(pf_incr_sm) <- rownames(falciparum_counts)

time2bins <- function(t, bins, width) { return(which( (bins-width) < t & (bins+width) > t)) }
pb_to_bins <- sapply(pseudotimes$pb, time2bins, bins=bins, width=width)
pk_to_bins <- sapply(pseudotimes$pk, time2bins, bins=bins, width=width)
pf_to_bins <- sapply(pseudotimes$pf, time2bins, bins=bins, width=width)

# Smoothed trans_rates
sm_fun_1d <- function(b, rate, time, width){
        subset <- time > b-width/2 & time < b+width/2;
        if (sum(subset) > 1) {
                return(mean(rate[subset]))
        } else if (sum(subset) == 1) {
                return(rate[subset])
        } else {
                return(NA);
        }
}
pb_rate_sm <- sapply(bins, sm_fun_1d, rate=trans_rates$pb, time=pseudotimes$pb, width=width)
pk_rate_sm <- sapply(bins, sm_fun_1d, rate=trans_rates$pk, time=pseudotimes$pk, width=width)
pf_rate_sm <- sapply(bins, sm_fun_1d, rate=trans_rates$pf, time=pseudotimes$pf, width=width)

smoothed_rates <- list(pb=pb_rate_sm, pk=pk_rate_sm, pf=pf_rate_sm)

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

pdf(paste("Pk_to_Pb_", dat_version, "_DTW_heatmap.pdf", sep=""), width=6, height=6)#, units="in", res=300)
heatmap.2(mat_bk, Rowv=FALSE, Colv=FALSE, trace="none", dendrogram="none", scale="column", col=colorRampPalette(c("grey90", "black"))(30), key.title="", key.xlab="Distance", xlab="P. bergei (t)", ylab="P. knowlesi (t)")
dev.off()

test_bk <- dtw(mat_bk, keep.internals=TRUE) 
test_bk$query <- pk_rate_sm
test_bk$reference <- pb_rate_sm

pdf(paste("Pk_to_Pb_", dat_version, "_DTW.pdf", sep=""), width=5.5, height=5.5)#, units="in", res=300)
my_dtwPlotTwoWay(test_bk, xlab="Time", ylab="Mean increase in expression", main="Pk vs Pb", yaxt="n")
dev.off()

pdf(paste("Pk_to_Pb_", dat_version, "_DTW_densityContour.pdf", sep=""), width=4, height=4)#, units="in", res=300)
my_dtwPlotDensity(test_bk, normalize=TRUE, xlab="Pk Time Index", ylab="Pb Time Index")
dev.off()

# index 1 = query index
# index 2 = reference index
# b vs f
mat_bf <- matrix(0, nrow=ncol(pf_sm_to_pb), ncol=ncol(pb_sm_to_pf));
for (i in 1:nrow(mat_bf)) {
	for(j in 1:ncol(mat_bf)) {
		mat_bf[i,j] <- d_bf(i,j)
	}
}

pdf(paste("Pf_to_Pb_", dat_version, "_DTW_heatmap.pdf", sep=""), width=6, height=6)#, units="in", res=300)
heatmap.2(mat_bf, Rowv=FALSE, Colv=FALSE, trace="none", dendrogram="none", scale="column", col=colorRampPalette(c("grey90", "black"))(30), key.title="", key.xlab="Distance", xlab="P. bergei (t)", ylab="P. falciparum (t)")
dev.off()

test_bf <- dtw(mat_bf, keep.internals=TRUE) 
test_bf$query <- pf_rate_sm
test_bf$reference <- pb_rate_sm
#deltaQs <- vector(length=length(unique(test_bf$index1)));
#names(deltaQs) <- colnames(pf_sm_to_pb)
#for (i in unique(test_bf$index1)) {
#	Refs <- unique(test_bf$index2[test_bf$index1==i])
#	deltaQ <- mean(bins[Refs]) - bins[i]
#	deltaQs[i] <- deltaQ
#}

#pf_shifts <- lapply(pf_to_bins, function(b, bins){mean(deltaQs[names(bins)[b]], na.rm=TRUE)}, bins=bins)

#pf_pseudotime_shifted <- pseudotimes$pf+unlist(pf_shifts)
# Just in case:
#pf_pseudotime_shifted[pf_pseudotime_shifted < 0] <- 0;
#pf_pseudotime_shifted[pf_pseudotime_shifted > 2*pi] <- 2*pi;

pdf(paste("Pf_to_Pb_", dat_version, "_DTW.pdf", sep=""), width=5.5, height=5.5)#, units="in", res=300)
my_dtwPlotTwoWay(test_bf, scalex=3.5,  xlab="Time", ylab="Mean increase in expression", main="Pf vs Pb")
dev.off()

pdf(paste("Pf_to_Pb_", dat_version, "_DTW_densityContour.pdf", sep=""), width=4, height=4)#, units="in", res=300)
my_dtwPlotDensity(test_bf, normalize=TRUE, xlab="Pf Time Index", ylab="Pb Time Index")
dev.off()

### Significance vs Permutations ###
set.seed(2310)
perm_ds <- vector()
start_time <- Sys.time()
for (r in 1:10000) {
	perm_row <- gtools::permute(1:nrow(mat_bk))
	perm_col <- gtools::permute(1:ncol(mat_bk))
	d <- dtw(mat_bk[perm_row, perm_col])$normalizedDistance
	perm_ds <- c(perm_ds, d)
}
end_time <- Sys.time()
saveRDS(perm_ds, paste("pb_vs_pk_permed_", dat_version, "_DTW_distance.rds", sep=""))

pdf(paste("Pk_to_Pb_", dat_version, "_DTW_PermDist.pdf", sep=""), width=4, height=4)#, units="in", res=300)
hist(c(perm_ds, test_bk$normalizedDistance), breaks=100, col="grey65", main="Pb-Pk vs Permutations", xlab="Distance")
arrows(test_bk$normalizedDistance, 100, test_bk$normalizedDistance, 0, col="red", lwd=3, length=0.15)
dev.off()

set.seed(4719)
perm_ds <- vector()
start_time <- Sys.time()
for (r in 1:10000) {
	perm_row <- gtools::permute(1:nrow(mat_bf))
	perm_col <- gtools::permute(1:ncol(mat_bf))
	d <- dtw(mat_bf[perm_row, perm_col])$normalizedDistance
	perm_ds <- c(perm_ds, d)
}
end_time <- Sys.time()
saveRDS(perm_ds, paste("pb_vs_pf_permed_", dat_version, "_DTW_distance.rds", sep=""))


pdf(paste("Pf_to_Pb_", dat_version, "_DTW_PermDist.pdf", sep=""), width=4, height=4)#, units="in", res=300)
hist(c(perm_ds, test_bf$normalizedDistance), breaks=100, col="grey65", main="Pb-Pf vs Permutations", xlab="Distance")
arrows(test_bf$normalizedDistance, 100, test_bf$normalizedDistance, 0, col="red", lwd=3, length=0.15)
dev.off()









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
#pdf("DTW_Three_Spp_Wave_Lines_w_Sig.pdf", width=7, height=9)
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
plot_piecewise(pseudotimes$pb, bergei_up_s, splits)

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
plot_piecewise(pseudotimes$pk, knowlesi_up_s, splits)

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
plot_piecewise(pseudotimes$pf, falciparum_up_s, splits)

#dev.off()

