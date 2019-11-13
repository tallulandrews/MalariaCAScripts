# My old fitting/residual method was to take the intersection of the top X% for x & y and force through 0,0. requiring that at least 10 cells be in said intersection.

# Read in various counts
exonic <- read.table("Exonic_Reads.out")
intronic <- read.table("Intronic_Reads.out")
spliced <- read.table("Spliced_Reads.out")

# Extract intron/exon info
extract_feature_data <- function(mat) {
	mat <- mat[-nrow(mat),]
	len <- unlist(mat[,1])
	mat <- mat[,-1]
	feat <- strsplit(rownames(mat), ":")
	feat <- t(sapply(feat, function(a) {c(a[1], a[length(a)-1], a[length(a)])}))
	#feat <- matrix(unlist(feat), ncol=4, byrow=TRUE)
	feat <- cbind(feat,len)
	colnames(feat) <- c("gene", "st", "end", "length")
	return(list(info=feat, MAT=mat))
}
exonic <- extract_feature_data(exonic)
exon.info <- exonic$info
exon.mat <- exonic$MAT
intronic <- extract_feature_data(intronic)
intron.info <- intronic$info
intron.mat <- intronic$MAT
spliced <- extract_feature_data(spliced)
spliced.mat <- spliced$MAT
spliced.info <- spliced$info

exon.mat <- exon.mat+spliced.mat

# Read in existing data
dat1 <- read.table("24829_1_1_intron_antisense.dat")
dat2 <- read.table("24829_1_2_intron_antisense.dat")
dat3 <- read.table("24829_1_3_intron_antisense.dat")

# Match up introns & exons
reorder <- order(exon.info[,1], as.numeric(exon.info[,2]));
exon.mat <- exon.mat[reorder,]
exon.info <- exon.info[reorder,]

reorder <- order(intron.info[,1], as.numeric(intron.info[,2]));
intron.mat <- intron.mat[reorder,]
intron.info <- intron.info[reorder,]

reorder <- order(spliced.info[,1], as.numeric(spliced.info[,2]));
spliced.mat <- spliced.mat[reorder,]
spliced.info <- spliced.info[reorder,]

reorder <- order(as.character(dat1[,1]), dat1[,2])
dat1 <- dat1[reorder,]
dat2 <- dat2[reorder,]
dat3 <- dat3[reorder,]

intron_keep <- intron.info[,1] %in% dat1[,1]
intron.matched.info <- intron.info[intron_keep,]
intron.matched.mat <- intron.mat[intron_keep,]

# Read in cell metadata
get_metadata <- function(rds_file) {
	require("scater")
	obj <- readRDS(rds_file)
	obj <- obj[,obj$use]
	colours <- as.character(obj$col)
	type <- as.character(obj$ShortenedLifeStage)
	if ("purple" %in% colours) {
		type[colours == "purple"] <- "male"
		type[colours == "purple4"] <- "female"
	}
	names(colours) <- colnames(obj);
	names(type) <- colnames(obj);
	return(list(cell_col=colours, cell_type=type));
}
Blood <- get_metadata("../allBLOOD.rds")	
IDC <- get_metadata("../allIDC.rds")	
SPZ <- get_metadata("../SPZ.rds")	
OOKOO <- get_metadata("../OOKOO.rds")	
EEF <- get_metadata("../EEF.rds")	

all_col <- c(Blood$cell_col, SPZ$cell_col, OOKOO$cell_col, EEF$cell_col)
all_type <- c(Blood$cell_type, SPZ$cell_type, OOKOO$cell_type, EEF$cell_type)

# Cell_name_matching
qc_cells1 <- read.table("../allBLOODgoodcells.csv", stringsAsFactors=FALSE, sep=",", header=T)
qc_cells2 <- read.table("../allIDCgoodcells.csv", stringsAsFactors=FALSE, sep=",", header=T)
qc_cells3 <- read.table("../SPZgoodcells.csv", stringsAsFactors=FALSE, sep=",", header=T)
qc_cells4 <- read.table("../OOKOOgoodcells.csv", stringsAsFactors=FALSE, sep=",", header=T)
qc_cells5 <- read.table("../EEFgoodcells.csv", stringsAsFactors=FALSE, sep=",", header=T)

qc_cells <- rbind(qc_cells1, qc_cells2, qc_cells3, qc_cells4, qc_cells5)
qc_cells$id <- paste("X", qc_cells[,2], "_", qc_cells[,3], sep="")
qc_cells <- unique(qc_cells);

do_cell_qc <- function(mat) {
	file_ids <- colnames(mat);
	file_ids <- sub("trim", "", file_ids)	
	file_ids <- sub("_sorted", "", file_ids)	
	keep_cells <- file_ids %in% qc_cells$id
	mat <- mat[, keep_cells]
	file_ids <- file_ids[keep_cells]
	cell_names <- qc_cells[match(file_ids,qc_cells$id), 1]
	colnames(mat) <- cell_names;
	return(mat);
}
exon.mat <- do_cell_qc(exon.mat)
intron.mat <- do_cell_qc(intron.mat)
intron.matched.mat <- do_cell_qc(intron.matched.mat)
spliced.mat <- do_cell_qc(spliced.mat)


# Normalization

SF <- colSums(exon.mat)+colSums(intron.mat)+colSums(spliced.mat);
overall_scale <- median(SF);

#exon.mat <- exon.mat/as.numeric(exon.info[,4])*100
#intron.mat <- intron.mat/as.numeric(intron.info[,4])*100
#intron.matched.mat <- intron.matched.mat/as.numeric(intron.matched.info[,4])*100
#spliced.mat <- spliced.mat/as.numeric(spliced.info[,4])*100

#SF <- colSums(exon.mat)+colSums(intron.mat)+colSums(spliced.mat);

exon.mat <- t( t(exon.mat)/SF*overall_scale)
intron.mat <- t( t(intron.mat)/SF*overall_scale)
intron.matched.mat <- t( t(intron.matched.mat)/SF*overall_scale)
spliced.mat <- t( t(spliced.mat)/SF*overall_scale)


# Aggregate by gene
require("CellTypeProfiles")
my_col_sum_aggregate <- function(mat, groups) {
	MAT <- as.matrix(mat)
	x <- split(seq(nrow(MAT)), groups)
	results <- t(sapply(x, function(a) if(length(a) > 1) {colSums(MAT[a,])} else {MAT[a,]}))
	return(results)
}

unspliced_by_gene <- my_col_sum_aggregate(intron.matched.mat, intron.matched.info[,1])
spliced_by_gene <- my_col_sum_aggregate(spliced.mat, spliced.info[,1])
exonic_by_gene <- my_col_sum_aggregate(exon.mat, exon.info[,1])

spliced_by_gene <- spliced_by_gene[rownames(spliced_by_gene) %in% rownames(unspliced_by_gene),]
exonic_by_gene <- exonic_by_gene[rownames(exonic_by_gene) %in% rownames(unspliced_by_gene),]

# Aggregate by type
all_type <- all_type[match(colnames(intron.mat), names(all_type))]
all_col <- all_col[match(colnames(intron.mat), names(all_col))]

intron_matched_by_type <- my_row_mean_aggregate(intron.matched.mat, all_type)
intron_by_type <- my_row_mean_aggregate(intron.mat, all_type)
exon_by_type <- my_row_mean_aggregate(exon.mat, all_type)
spliced_by_type <- my_row_mean_aggregate(spliced.mat, all_type)

# Get antisense intron
get_antisense <- function(g, dat) {
	# First highest ratio of anti:sense
	# Second highest abs(anti)
	g_rows <- which(dat[,1] %in% g & dat[,4] >=5);
	if (length(g_rows) == 0) {return(NA)}
	if (length(g_rows) == 1) {
		return(g_rows)
	} else {
		this <- dat[g_rows,]
		ratio <- this[,4]/this[,3]
		anti <- which(ratio == max(ratio));
		if (length(anti)==1) {return(g_rows[anti])}
		anti <- which(this[anti,4] == max(this[anti,4]));
		return(g_rows[anti[1]]);
	}	
}

### For Anti-sense Require: Correlation across types < -0.2
get_cors_type <- function(g, dat) {
	E<-exon_by_type[exon.info[,1] == g,]
	if (nrow(E) > 1) {
		E <- colSums(E);
	}
	anti_intron <- get_antisense(g, dat)
	if (!is.na(anti_intron)) {
		Anti <- intron_matched_by_type[anti_intron,]
	} else {
		Anti <- rep(0, times=ncol(intron_matched_by_type))
	}
	cor(Anti, E)
}

get_cors <- function(g, dat) {
	E<-exon.mat[exon.info[,1] == g,]
	if (nrow(E) > 1) {
		E <- colSums(E);
	}
	anti_intron <- get_antisense(g, dat)
	if (!is.na(anti_intron)) {
		Anti <- intron.matched.mat[anti_intron,]
	} else {
		Anti <- rep(0, times=ncol(intron.matched.mat))
	}
	cor(Anti, E)
}

alldat <- dat1;
alldat[,3] <- dat1[,3]+dat2[,3]+dat3[,3]
alldat[,4] <- dat1[,4]+dat2[,4]+dat3[,4]

E_Anti_cors <- sapply(unique(intron.matched.info[,1]), get_cors_type, alldat)
E_Anti_cors[is.na(E_Anti_cors)] <- 0;

potential <- names(E_Anti_cors)[E_Anti_cors < -0.2]

### For Antisense require: Most intronic reads come from anti-intron
anti_expression <- function(g, dat) {
	anti_intron <- get_antisense(g, dat)
	if (!is.na(anti_intron)) {
		return(intron.matched.mat[anti_intron,])
	} else {
		return(rep(0, times=ncol(intron.matched.mat)));
	}
}

is.anti <- function(g, dat) {
	E<-exon.mat[exon.info[,1] == g,]
	if (sum(exon.info[,1] == g) > 1) {
		E <- colSums(E);
	}
	anti_intron <- get_antisense(g, dat)
	if (!is.na(anti_intron)) {
		Anti <- intron.matched.mat[anti_intron,]
	} else {
		#Anti <- rep(0, times=ncol(intron.matched.mat))
		return(rep(FALSE, times=ncol(intron.matched.mat)));
	}
	I <- intron.mat[intron.info[,1]==g,]
	if (sum(intron.info[,1]==g) > 1) {
		I <- colSums(I);
	} else {
		return(rep(FALSE, times=ncol(intron.matched.mat)));
	}
	ratio <- unlist(Anti)/(unlist(I));
	ratio[is.na(ratio)] <- 0;
	threshold <- mean(ratio[I > 0])
	return(ratio >= threshold);
}

isAnti_by_cell <- t(sapply(potential, is.anti, alldat))
isAnti_by_type <- my_row_mean_aggregate(isAnti_by_cell, all_type)
Anti_by_cell <- t(sapply(potential, anti_expression, alldat))


write.table(isAnti_by_type, file=paste(prefix, "prop_AntiSense_by_type.csv", sep="_"))


