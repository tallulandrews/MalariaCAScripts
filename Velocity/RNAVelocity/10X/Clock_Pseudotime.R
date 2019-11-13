# Arguments:
# pcs = matrix of principle components, rows = cell, columns = PCs
# start.cell = index/row of the starting cell.
# direction = does pseudotime go clockwise or counter-clockwise around the ring
# suppress.plot = plot the fitted ellipse and pseudotime.
# include.radial.distance = account for the distance from the centre of the 
#	ellipse when calculating pseudotime, not just the angle.
# window.size = number of cells in a window for smoothing the radial distances.
# smoothing.factor = parameter for amount of smoothing of the radial distances (high = more)

# Output:
# fit.pseudotime = estimated pseudotime
# fit.ellipse = ellipse parameters
# fit.x = fitted x-coordinate of the cell in pca space (radial distance only)
# fit.y = fitted y-coordinate of the cell in pca space (radial distance only)
# fit.angles = fitted angles for each cell (if no radial distance this is the same as fit.pseudotime.

clock_pseudotime <- function(pcs, start.cell=1, direction=c("clockwise", "counter"), suppress.plot=FALSE, include.radial.distance=FALSE, window.size=20, smoothing.factor=0.7) {

        if (length(start.cell) != 1) {stop("Error: Please supply exactly one starting cell.")}
        if (class(start.cell) != "numeric") {stop("Error: Supply starting cell as a numeric index.")}
        if (ncol(pcs) != 2) {warning("Warning: currently only 2D ring fitting is implemented. Using first two columns of pcs.")}

	 require("conicfit")
        # Fit Centre
        fit <- EllipseDirectFit(pcs)
        fit <- AtoG(fit)[[1]]
        ellipse.centre <- c(fit[1,1], fit[2,1])

        # Calculate Angles
        angles <- atan2(pcs[,2]-ellipse.centre[2], pcs[,1]-ellipse.centre[1]);
	angles.shifted <- angles-angles[start.cell]
	angles.shifted[angles.shifted < 0] <- 2*pi-angles[start.cell] +
                angles.shifted[angles.shifted < 0] + angles[start.cell]

        if (direction[1] == "counter") {
		angles.shifted <- 2*pi-angles.shifted.
        } 
	fit.pseudotime <- angles.shifted
	fit.x <- NULL
	fit.y <- NULL
	cells.ordered <- c(1:length(angles.shifted))[order(angles.shifted)]

	if (include.radial.distance) {
		# Calculate smoothed distance to centre
        	cell.angles <- angles.shifted
        	smoothed.radius <- rep(0, times=length(cell.angles))
        	fit.pseudotime <- rep(0, times=length(cell.angles))
        	fit.x <- rep(0, times=length(cell.angles))
        	fit.y <- rep(0, times=length(cell.angles))
        	flank <- ceiling(window_size/2)
        	for(i in 1:length(angles.shifted)) {
        	        pos <- cells.ordered[i]
        	        if (pos-flank < 1) {
        	                window <- c( (length(cells.ordered)-(flank-pos)):length(cells.ordered), 1:pos, (pos+1):(pos+flank))
                	} else if (pos+flank > length(cells.ordered)) {
                	        window <- c( (pos-flank):(pos-1), pos:length(cells.ordered),
                                        1:( flank-(length(cells.ordered) - pos) ) )
                	} else {
                	        window<- c( (pos-flank):(pos+flank) )
                	}

                	window.cells <- cells.ordered[window]
                	smoothed.radius[pos] <- median( sqrt((pcs[window.cells,1]-ellipse.centre[1])^2 + (pcs[window.cells,2]-ellipse.centre[2])^2  ))
        	}

        	smoothed.radius[cells.ordered] <- smooth.spline(as.vector(smoothed.radius[cells.ordered]), spar=smoothing.factor)$y
		# account for distance to centre in pseudotime
	        theta <- atan2(pcs[,2]-ellipse.centre[2], pcs[,1]-ellipse.centre[1])
	        fit.x <- smoothed.radius*cos(theta) + ellipse.centre[1]
	        fit.y <- smoothed.radius*sin(theta) + ellipse.centre[2]
		for(i in 1:length(angles.shifted)) {
                	pos <- cells.ordered[i]
                	if (i == 1) {
                	        fit.pseudotime[pos] <- 0
                	} else {
                	        pos_prev <- cells.ordered[i-1];
                	        fit.pseudotime[pos] <- fit.pseudotime[pos_prev]+
                                	dist(rbind(c(fit.x[pos],fit.y[pos]),
                                        c(fit.x[pos_prev],fit.y[pos_prev])));
                	}
        	}
	}

	if (!suppress.plot) {
                require("RColorBrewer")
                vals <- cut(fit.pseudotime, breaks=7)
                cols <- brewer.pal(8, "Greys")[2:8]
                plot(pcs[,1], pcs[,2], col=cols[vals], pch=16)
                plot_fit <- calculateEllipse(fit[1,1], fit[2,1], fit[3,1], fit[4,1], fit[5,1], 300)
                lines(plot_fit, lty=3, lwd=1.5, col="grey25")
                points(pcs[start.cell,1], pcs[start.cell,2], pch=1, col="red")
		if (include.radial.distance) {
                	lines(fit.x[cells.ordered], fit.y[cells.ordered], 
				col="blue", lty=1, lwd=2.5)
		}
        }


	return(list(pseudotime=fit.pseudotime, fit.ellipse=fit, fit.x=fit.x, fit.y=fit.y, fit.angle=angle.shifted))
}





