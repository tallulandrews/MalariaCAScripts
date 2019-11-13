# Examining transcriptional waves using RNA velocity through IDC

## Smartseq2
require("scater")
dat_SS <- readRDS("../NewBatches/IDC_RNAVelo_Counts.rds")
velo_SS <- readRDS("../NewBatches/IDC_RNAVelo_Estimates.rds")
ann_SS <- readRDS("../NewBatches/allIDC.rds")
pca_SS <- readRDS("../NewBatches/IDC_RNAVelo_PCAobj.rds")

## 10X velocyto
dat_10X <- readRDS("bergei_10X_velodata.rds")
velo_10X <- readRDS("bergei_10X_veloEst.rds")

# Testing
pcs <- pca_SS$epc@scores[,1:2]
start_cell=which(ann_SS$ShortenedLifeStage =="Ring")

# fit an ellipse
require("conicfit")
fita <- EllipseDirectFit(dat_10X$ann[,c("pc1", "pc2")])
fitg <- AtoG(fita)[[1]]
plot_fit <- calculateEllipse(fitg[1,1], fitg[2,1], fitg[3,1], fitg[4,1], fitg[5,1], 300)

plot(dat_10X$ann[,c("pc1", "pc2")])
points(plot_fit, col="red")

res <-ResidualsG(as.matrix(dat_10X$ann[,c("pc1", "pc2")]), fitg)

# center the ellipse
tmp_xes <- dat_10X$ann[,"pc1"]-fitg[1,1]
tmp_yes <- dat_10X$ann[,"pc2"]-fitg[2,1]
plot(tmp_xes, tmp_yes)
# Calculate angle of each cell.
angles <- atan2(tmp_yes,tmp_xes)
tmp_col <- brewer.pal(8,"Greys")[2:8]
angles <- pi-angles
plot(tmp_xes, tmp_yes, col=tmp_col[as.numeric(cut(angles, breaks=7))], pch=16)

start_cell <- which(round(tmp_xes)==-3 & round(tmp_yes)==-6 & dat_10X$ann$cell_name=="ring")
start_cell <- which(angles == min(angles[start_cell]))
points(tmp_xes[start_cell], tmp_yes[start_cell], col="red", pch=16)

converted <- angles-angles[start_cell]
converted[converted < 0] <- max(converted)+abs(min(converted[converted < 0]))+converted[converted < 0]


fit_Ring <- function(pcs, start_cell=1, direction=c("clockwise", "counter"), suppress.plot=FALSE, window_size=20, smoothing.factor=0.7) {

	if (length(start_cell) != 1) {stop("Error: Please supply exactly one starting cell.")}
	if (ncol(pcs) != 2) {warning("Warning: currently only 2D ring fitting is implemented. Using first two columns of pcs.")}

	require("conicfit")
	# Fit Centre
	fit <- EllipseDirectFit(pcs)
	fit <- AtoG(fit)[[1]]
	fit_centre <- c(fit[1,1], fit[2,1])

	# Calculate Angles
	angles <- atan2(pcs[,2]-fit_centre[2], pcs[,1]-fit_centre[1]);
	if (direction[1] == "clockwise") {
		angles <- pi-angles;
	} else {
		angles <- angles+pi
	}
	converted <- angles-angles[start_cell]
	converted[converted < 0] <- max(converted)+
		abs(min(converted[converted < 0]))+
		converted[converted < 0]

	# Calculate smoothed distance to centre
	ordered_pts <- c(1:length(converted))[order(converted)]
	fit_angles <- converted
	fit_radius <- rep(0, times=length(fit_angles))
	fit_pseudotime <- rep(0, times=length(fit_angles))
	fit_x <- rep(0, times=length(fit_angles))
	fit_y <- rep(0, times=length(fit_angles))
	flank <- ceiling(window_size/2)
	for(i in 1:length(converted)) {
		pos <- ordered_pts[i]
		if (pos-flank < 1) {
			window <- c( (length(ordered_pts)-(flank-pos)):length(ordered_pts),
					 1:pos, (pos+1):(pos+flank))
		} else if (pos+flank > length(ordered_pts)) {
			window <- c( (pos-flank):(pos-1), pos:length(ordered_pts),
					1:( flank-(length(ordered_pts) - pos) ) )
		} else {
			window<- c( (pos-flank):(pos+flank) )
		}

		window_pts <- ordered_pts[window]
		fit_radius[pos] <- median( sqrt((pcs[window_pts,1]-fit_centre[1])^2 + (pcs[window_pts,2]-fit_centre[2])^2  ))
	}

	fit_radius[ordered_pts] <- smooth.spline(as.vector(fit_radius[ordered_pts]), spar=smoothing.factor)$y
	theta <- atan2(pcs[,2]-fit_centre[2], pcs[,1]-fit_centre[1])
	fit_x <- fit_radius*cos(theta) + fit_centre[1]
	fit_y <- fit_radius*sin(theta) + fit_centre[2]

	for(i in 1:length(converted)) {
		pos <- ordered_pts[i]
		if (i == 1) {
			fit_pseudotime[pos] <- 0
		} else {
			pos_prev <- ordered_pts[i-1];
			fit_pseudotime[pos] <- fit_pseudotime[pos_prev]+
				dist(rbind(c(fit_x[pos],fit_y[pos]), 
					c(fit_x[pos_prev],fit_y[pos_prev])));
		}
	}
	
	if (!suppress.plot) {
		require("RColorBrewer")
		vals <- cut(fit_pseudotime, breaks=7)
		cols <- brewer.pal(8, "Greys")[2:8]
		plot(pcs[,1], pcs[,2], col=cols[vals], pch=16)
		plot_fit <- calculateEllipse(fit[1,1], fit[2,1], fit[3,1], fit[4,1], fit[5,1], 300)
		lines(plot_fit, lty=3, lwd=1.5, col="grey25")
		points(pcs[start_cell,1], pcs[start_cell,2], pch=1, col="red")
		lines(fit_x[ordered_pts], fit_y[ordered_pts], col="blue", lty=1, lwd=2.5)
	}
	return(list(pseudotime=fit_pseudotime, fit_x=fit_x, fit_y=fit_y, pseudo_order=ordered_pts, orig_ellipse=fit, orig_angle=angle, rot_angle=converted))


}
