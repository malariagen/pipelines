# Load a library that will allow padding 0s in front of a number
library(stringr)
library(parallel)

# This is the minimum number of consecutive duplicated windows required to call a duplication
min.dup.length <- 5
# This is the threshold variance in coverage above which we remove a sample
threshold.variance <- 0.2
# This is the coverage window that was used to calculate the coverage for the samples we are loading. 
# This parameter should not be changed unless using input data calculated with a different window size.
coverage.window <- 300 
# Set the minimum proportion of windows in a gene that need to have some coverage data before making a duplication call.
min.win.cov <- 0.5
# Set the likelihood ratio threshold at which we will accept a CNV
lik.thresh <- 10000

arg.values <- commandArgs(trailingOnly=T)
cat('Arguments:', arg.values, '\n', sep = '\n')

# Get the chrosome
acceptable.chromosomes <- c('2L', '2R', '3L', '3R', 'X')
chrom <- arg.values[1]
if (!(chrom %in% acceptable.chromosomes)){
	stop(paste('Chromosome should be one of: ', paste(acceptable.chromosomes, collapse = ', ')))
}

# Get the names of the samples to analyse
sample.list <- arg.values[2]
sample.names <- setNames(nm = read.table(sample.list, stringsAsFactors = F)[,1])

# Load the coverage variance data
cov.var.file <- arg.values[3]
cov.var <- read.table(cov.var.file, header = T, sep = '\t', row.names = 1)
# Identify samples with high coverage variance for exclusion
high.var.samples <- intersect(sample.names, rownames(cov.var)[cov.var$autosomes > threshold.variance])
good.var.samples <- setdiff(sample.names, rownames(cov.var)[cov.var$autosomes > threshold.variance])

# Load up the information on gene and exon positions. 
gene.pos.file <- arg.values[4]
gene.positions <- read.table(gene.pos.file)
gene.positions <- gene.positions[gene.positions$Chrom == chrom, ]

# Get the detox genes file
detox.genes.file <- arg.values[5]
detox.genes <- read.table(detox.genes.file, stringsAsFactors = F)[,1]

# Get the working directory
working.directory <- arg.values[6]

# Get the number of cores
n.cores <- arg.values[7]

# Get the output folder
output.folder <- arg.values[8]
dir.create(output.folder, showWarnings = FALSE)

# If this is the sex chromosome, we also need the metadata to determine sex. 
	# These are the minimum copy number states required to consider a duplication (2 = "normal" diploid state).
{if (chrom == 'X'){
	metafile <- arg.values[9]
	meta <- read.table(metafile, sep = '\t', header = T, row.names = 1, quote = '', comment.char = '')[sample.names, ]
	threshold.copy.number <- c(3, 2)[(meta$sex_call == 'M') + 1]
}
else {
	threshold.copy.number <- 3
}}


# Function to identify runs of increased copy number state
find.full.dup.sequences <- function(cnv, threshold, n){
	# If cnv is a dataframe, extract the cnv vector
	if (is.data.frame(cnv))
		cnv.vector <- cnv$CNV
	else
		cnv.vector <- cnv
	current.run <- 0
	output <- matrix(0,0,2)
	for (i in 1:length(cnv.vector)){
		if (cnv.vector[i] >= threshold) {
			if (current.run == 0)
				start.pos <- i
			current.run <- current.run + 1
		}
		else{
			if (current.run >= n)
				output <- rbind(output, c(start.pos, i-1))
			current.run <- 0
		}
	}
	if (current.run >= n)
		output <- rbind(output, c(start.pos, i-1))
	return(output)
}


# Function to calculate the likelihood ratio of a duplication vs no duplication for a given CNV
# call
lik.ratio <- function(hmm, sample.variance, base.variance = 0.01, trans.prob = 0.00001, ploidy = 2){
	# Get the probability of not changing state and the ratio of the two probabilities
	non.trans.prob <- 1 - 12*trans.prob
	# Get the number of transitions in the cnv states
	transitions <- sum((hmm$CNV[-1] - hmm$CNV[-nrow(hmm)]) != 0)
	# get the ratio of transition probabilities
	trans.ratio <- (trans.prob^transitions)/(non.trans.prob^transitions)
	# get the likelihood of the null model. 
	lik.null <- dnorm(hmm$Normalised_coverage, ploidy, sqrt(base.variance + (sample.variance/2)*ploidy))
	# get the likelihood of the duplication model. 
	lik.dup <- dnorm(hmm$Normalised_coverage, hmm$CNV, sqrt(base.variance + (sample.variance/2)*hmm$CNV))
	# If the likelihood of the duplicated model at any given window is 0, then coverage must have been so high
	# that even CNV state 12 could not capture it, so we report a very large likelihood ratio
	if (any(lik.dup == 0))
		return(Inf)
	# Otherwise, get the product of the ratios multiplied by the transition probability
	else
		return(trans.ratio * prod(lik.dup/lik.null))
}

# Function to identify the genetic features encompass in a CNV
find.features.encompassed <- function(pos.table, genepos.table, filtered.windows, cov.win){
	output.table <- data.frame(matrix(0,0,3, dimnames = list(character(),c('start', 'end', 'Features'))))
	if (nrow(pos.table) == 0){
		return(output.table)
	}
	for (i in 1:nrow(pos.table)){
		# Find the genes whose post-filtering windows are all included in this duplication
		gene.req <- genepos.table[,'start'] >= pos.table[i,1] & genepos.table[,'end'] <= pos.table[i,2]
		# Gene that do not have enough covered windows to be considered will produce NAs in the above line. Turn 
		# these to FALSE
		gene.req[is.na(gene.req)] <- F
		# The regions where both of these requirements are met represent a gene encompassed in the duplication
		genes.encompassed <- genepos.table[gene.req,]
		# If there is at least one overlap, add the pos.table region to the output table and record the name of the 
		# genepos.table regions that overlapped it
		if (nrow(genes.encompassed)){
			region.names <- paste(rownames(genes.encompassed), collapse = ";")
			new.row <- data.frame(start = pos.table[i,1], end = pos.table[i,2], Features = region.names, stringsAsFactors = F)[1,]
			if (!is.null(rownames(pos.table)))
				output.table[rownames(pos.table)[i],] <- new.row
			else
				output.table <- rbind(output.table, new.row)
		}
	}
	# Now convert the indices to genomic positions in the output table
	output.table[,1] <- filtered.windows[output.table[,1]]
	output.table[,2] <- filtered.windows[output.table[,2]] + cov.win
	return(output.table)
}

# Function that will go through all of the duplications in a list and merge them if they have similar start 
# and end positions. Let's write it such that the function can be run on either indices or positions, making 
# the permissiveness around the start and end points an argument, so that appropriate values can be given 
# depending on whether indices or positions are used. 
merge.duplications <- function(dup.list, leeway){
	# Turn the matrices that make up the input list into dataframes
	dup.list.out <- lapply(dup.list, function(x) data.frame(start = x[,1], end = x[,2], CNV_code = 0))
	# Here is the code we will assign the next CNV group we find
	current.CNV.code = 1
	# Now go through every dup, looking for matching ranges
	for (i in 1:length(dup.list.out)){
		cat('Analysing sample ', names(dup.list.out)[i], '.\n', sep = '')
		for (k in 1:nrow(dup.list.out[[i]])){
			this.start <- dup.list.out[[i]][k,1]
			this.end <- dup.list.out[[i]][k,2]
			# Go through every previous sample, finding matching duplications. We record matching duplications 
			# by the sample and row in which they are found, and their current CNV_code.
			matching.dups <- matrix(0,0,3)
			for (j in (1:i)[-i]){
				# find any dups that match the start and end points of this one, including the leeway
				match.start <- (dup.list.out[[j]][,1] >= this.start - leeway) & (dup.list.out[[j]][,1] <= this.start + leeway)
				match.end <- (dup.list.out[[j]][,2] >= this.end - leeway) & (dup.list.out[[j]][,2] <= this.end + leeway)
				matching.dup.index <- which(match.start & match.end)
				# We should never find more than one matching dup in a given sample. Check that this is true
				if (length(matching.dup.index) > 1) 
					stop(paste('Fail. Found more than one matching duplications for i = ', i, ', k = ', k, ', j = ', j, '.'))
				else if (length(matching.dup.index) == 1)
					matching.dups <- rbind(matching.dups, c(j, matching.dup.index, dup.list.out[[j]][matching.dup.index, 3]))
			}
			# Now we have identified all of the dups that match this one, we need to merge them all together
			# If there are no matching dups, then give the current CNV the next available code
			if (nrow(matching.dups) == 0){
				dup.list.out[[i]][k,3] <- current.CNV.code
				current.CNV.code <- current.CNV.code + 1
			}
			else {
				# If there is only one matching dup, just give the same code to the current dup
				matching.dup.codes <- sort(unique(matching.dups[,3]))
				if (length(matching.dup.codes) == 1)
					dup.list.out[[i]][k,3] <- matching.dups[1,3]
				# Otherwise, find the smallest of the matching codes and give all matching dups that code. We
				# also need to find all the dups that have the matching code (not just the ones that directly
				# matched this dup) and change the code for those dups. 
				else {
					smallest.code <- matching.dup.codes[1]
					cat('\tMerging CNV codes ', paste(matching.dup.codes, collapse = ', '), '.\n', sep = '')
					# Give this code to the current dup
					dup.list.out[[i]][k,3] <- smallest.code
					# For each other dup code, find all the samples that have that code and change it to the 
					# smallest code
					for (j in (1:i)[-i]){
						dup.list.out[[j]][dup.list.out[[j]][,3] %in% matching.dup.codes[-1], 3] <- smallest.code
					}
				}
			}
		}
	}
	dup.list.out
}

# Function to load up each file and store the results in a list
load.hmm.file <- function(sample.name){
	all.files <- list.files()
	# Find the file corresponding to this sample
	this.file <- grep(sample.name, all.files, value = T)
	if (length(this.file) == 0)
		stop(paste('No matching files were found for sample', sample.name))
	if (length(this.file) > 1)
		stop(paste('More than one matching file were found for sample', sample.name))
	# Load the file 
	cat('\t', this.file, '\n', sep = '')
	read.table(this.file, header = T)
}

# Function to filter samples by likelihood ratio
filter.by.likelihood.ratio <- function(cnv.table, hmm.table, sample.name, likelihood.threshold){
	filtered.table <- matrix(0,0,2)
	keeprow <- function(x)
		lik.ratio(hmm.table[x[1]:x[2], ], cov.var[sample.name, 'autosomes']) > likelihood.threshold
	filtered.table <- cnv.table[apply(cnv.table, 1, keeprow), ]
	filtered.table
}

# Function to load the data for a sample and detect raw CNVs within it
process.sample <- function(sample.name, threshold, n, likelihood.threshold){
	these.hmm <- load.hmm.file(sample.name)
	duplications.indices <- find.full.dup.sequences(these.hmm, threshold, n)
	goodlr.duplications.indices <- filter.by.likelihood.ratio(duplications.indices, these.hmm, sample.name, likelihood.threshold)
	goodlr.duplications.indices
}

# Set the working directory
setwd(working.directory)

cat('Processing HMM files\n')
# We'll use the first sample to get the list of posistions that remained after window filtering
hmm.positions <- load.hmm.file(sample.names[1])$Position
# Now load up all the samples and look for CNVs
goodlr.duplications.list.indices <- mcmapply(process.sample, sample.names, threshold = threshold.copy.number, n = min.dup.length, likelihood.threshold = lik.thresh, mc.cores = n.cores, SIMPLIFY = F)
# Because I am not 100% sure that mcmapply will output the data in the right order, the next line ensures the order
# is correct
goodlr.duplications.list.indices <- goodlr.duplications.list.indices[sample.names]

## Now create a single table listing all of these raw CNVs in all samples, including the high variance samples
num.raw.cnvs <- sapply(goodlr.duplications.list.indices, nrow)
full.raw.cnv.indices <- do.call(rbind, goodlr.duplications.list.indices)
full.raw.cnv.table <- data.frame(Sample = rep(names(num.raw.cnvs), num.raw.cnvs),
                                 Chrom = rep(chrom, sum(num.raw.cnvs)),
                                 # Get the CNV coordinates as positions rather than indices
                                 matrix(hmm.positions[full.raw.cnv.indices] + rep(c(0,coverage.window), each = nrow(full.raw.cnv.indices)), ncol = 2)
                                )
colnames(full.raw.cnv.table)[3:4] <- c('Start', 'Stop')
write.table(full.raw.cnv.table, paste(output.folder, '/full_raw_CNV_table_', chrom, '.csv', sep = ''), sep = '\t', row.names = F, quote = F)

# Now remove the samples with high variance
goodvar.goodlr.duplications.list.indices <- goodlr.duplications.list.indices[good.var.samples]

# Merge CNV clusters, but only use the samples with good variance values for this
cat('Merging CNV clusters.\n')
clustered.list <- merge.duplications(goodvar.goodlr.duplications.list.indices, 1)
save.image(paste('temp.Rdata', sep = ''))

# Order these data by CNV cluster, instead of by sample
combined.clustered.list <- do.call(rbind, clustered.list)
all.CNVs <- unique(combined.clustered.list[,3])
# The names of the all.CNVs vector should be the names that we will give to the CNVs. We pad them with at 
# least 4 characters, or more if required. 
names(all.CNVs) <- paste('CNV_', chrom, str_pad(all.CNVs, max(c(4, max(nchar(all.CNVs)))), pad = '0'), sep = '')
# List the data by CNV, instead of by sample.
get.cnv.cluster <- function(cnv, cnv.table){
	this.CNV.table <- subset(cnv.table, CNV_code == cnv)
	# Unlike phase2, some CNV codes end up existig twice in the same sample (when the range in start and 
	# end positions ends up large enough. So the sample names can no longer be the row names of the CNV 
	# table
	this.CNV.table$sample.name <- sub('\\..+', '', rownames(this.CNV.table))
	this.CNV.table[,c('sample.name', 'start', 'end')]
}
duplications.by.cluster <- lapply(all.CNVs, get.cnv.cluster, combined.clustered.list)
rm(all.CNVs)
rm(combined.clustered.list)

# Create a table showing the start and end points of these CNVs. For each start point, we take the median of 
# the start points of all the occurences of that CNV. We also want a measure of the uncertainty of the values
# of those medians, which we calculate as the range which each of the start and end points take.
summarise.CNVs <- function(CNV.cluster){
	start.point = ceiling(median(CNV.cluster$start))
	end.point = floor(median(CNV.cluster$end))
	start.05.quantile = floor(quantile(CNV.cluster$start, 0.05))
	start.95.quantile = ceiling(quantile(CNV.cluster$start, 0.95))
	end.05.quantile = floor(quantile(CNV.cluster$end, 0.05))
	end.95.quantile = ceiling(quantile(CNV.cluster$end, 0.95))
	num_samples = length(unique(CNV.cluster$sample.name))
	setNames(c(start.point, end.point, start.05.quantile, start.95.quantile, end.05.quantile, end.95.quantile, num_samples),
	         c('start', 'end', 'start.05.quantile', 'start.95.quantile', 'end.05.quantile', 'end.95.quantile', 'num_samples'))
}
CNV.table <- data.frame(t(sapply(duplications.by.cluster, summarise.CNVs)))

# For each duplications table, identify any genes that are encompassed by any duplications, with a minimum
# number of windows providing useable coverage data. 
# First, convert the positions of the gene table into indices of windows that passed filtering. 
determine.eligible.gene.positions <- function(this.gene.coords, positions, cov.win, min.cov){
	this.gene.start <- this.gene.coords[1]
	this.gene.end <- this.gene.coords[2]
	# Get the windows included in this gene
	this.gene.covered.windows <- which((positions >= this.gene.start) & (positions <= this.gene.end - cov.win))
	# Get the number of windows that would be included in this gene if none had been filtered
	this.gene.full.windows <- floor(this.gene.end/cov.win) - ceiling(this.gene.start/cov.win)
	# If there are no potential windows fully encompassed in the gene, declare the gene ineligible
	if (this.gene.full.windows <= 0)
		return(c(length(this.gene.covered.windows),
		         0,
		         NA, 
		         NA))
	else if (length(this.gene.covered.windows) >= (min.cov * this.gene.full.windows))
		return(c(length(this.gene.covered.windows),
		         this.gene.full.windows,
		         min(this.gene.covered.windows), 
		         max(this.gene.covered.windows)))
	else 
		return(c(length(this.gene.covered.windows),
		         this.gene.full.windows,
		         NA, 
		         NA))
}
eligible.gene.position.indices <- cbind(gene.positions[, 'Chrom', drop = F], data.frame(t(matrix(apply(gene.positions[, c('start', 'end')], 1, determine.eligible.gene.positions, hmm.positions, coverage.window, min.win.cov), nrow = 4, dimnames = list(c('covered.windows', 'full.windows', 'start', 'end'))))))

# Now look for encompassed genes
cat('Finding genes encompassed within CNVs\n')
# Note that the eligible.gene.position.indices table records the indices of the windows that are fully contained
# inside each gene. Thus, implicitly, we are only requiring that these windows are in a cnv in order to call a
# gene as duplicated. If the windows either side, overlapping the ends of the gene) are not in the cnv, this will
# not prevent the duplication call. 
CNV.encompass.table <- find.features.encompassed(CNV.table, genepos.table = eligible.gene.position.indices, filtered.windows = hmm.positions, cov.win = coverage.window)
# Create a list that captures the genes present in each CNV in a more effective way.
CNV.encompass.list <- strsplit(CNV.encompass.table$Features, ';')
names(CNV.encompass.list) <- rownames(CNV.encompass.table)
CNV.encompass.table$has.detox <- sapply(CNV.encompass.list, function(x) any(x %in% detox.genes))
# Create the full.CNV.table. Num_samples is not necessarily gambiae + arabiensis because there is one intermediate
# sample.
full.CNV.table <- data.frame(Chrom = chrom,
                             Start = hmm.positions[CNV.table$start],
                             End = hmm.positions[CNV.table$end] + coverage.window,
                             Start_90_quantile_range = hmm.positions[CNV.table$start.95.quantile] - hmm.positions[CNV.table$start.05.quantile],
                             End_90_quantile_range = hmm.positions[CNV.table$end.95.quantile] - hmm.positions[CNV.table$end.05.quantile],
                             Num_samples = CNV.table$num_samples,
                             Genes = CNV.encompass.table[rownames(CNV.table), 'Features'],
                             Has_detox = CNV.encompass.table[rownames(CNV.table), 'has.detox'],
							 row.names = rownames(CNV.table),
                             stringsAsFactors = F)
full.CNV.table$Genes[is.na(full.CNV.table$Genes)] <- ''
full.CNV.table$Has_detox[is.na(full.CNV.table$Has_detox)] <- F
for (s in good.var.samples)
	full.CNV.table[[s]] <- sapply(duplications.by.cluster, function(x) s %in% x$sample.name)
rm(CNV.encompass.table)

write.table(full.CNV.table, paste(output.folder, '/full_coverage_CNV_table_', chrom, '.csv', sep = ''), sep = '\t', col.names = NA, quote = F)

save.image(paste(output.folder, '/CNV_analysis_', chrom, '.Rdata', sep = ''))

