library(stringr)
library(stringi)
library(Biostrings)
library(parallel)

arg.values <- commandArgs(trailingOnly=T)
cat('Arguments:', arg.values, '\n', sep = '\n')

threshold.variance <- 0.2

sample.list <- arg.values[1]
sample.names <- setNames(nm = read.table(sample.list, stringsAsFactors = F)[[1]])

gene.regions.file <- arg.values[2]
all.gene.coordinates <- read.table(gene.regions.file, sep = '\t', header = T, row.names = 1)

meta.file <- arg.values[3]
meta <- read.table(meta.file, sep = '\t', header = T, row.names = 1, quote = '', comment.char = '')[sample.names, ]
expected.copy.number.on.X <- c(2, 1)[(meta$sex_call == 'M') + 1]

# Get the species calls
species.calls.file <- arg.values[4]
species.calls <- read.table(species.calls.file, sep = '\t', header = T, row.names = 1)

cov.var.file <- arg.values[5]
cov.var <- read.table(cov.var.file, header = T, sep = '\t', row.names = 1)
high.var.samples <- intersect(sample.names, rownames(cov.var)[cov.var$autosomes > threshold.variance])

coverage.folder <- arg.values[6]
diagnostic.reads.folder <- arg.values[7]
plotting.functions.file <- arg.values[8]
num.cores <- as.numeric(arg.values[9])

# Function to load up each file and store the results in a list
load.hmm.file <- function(sample.name, folder = '.'){
	all.files <- list.files(folder)
	# Find the file corresponding to this sample
	this.file <- paste(folder, grep(sample.name, all.files, value = T), sep = '/')
	if (length(this.file) == 0)
		stop(paste('No matching files were found for sample', sample.name))
	if (length(this.file) > 1)
		stop(paste('More than one matching file were found for sample', sample.name))
	# Load the file 
	cat('\t', this.file, '\n', sep = '')
	read.table(this.file, header = T)
}

# Function that will take a table of HMM output and shrink it to contain only a desired region
shrink.compact.hmm <- function(region.coords, input.hmm){
	subset(input.hmm, (Position >= region.coords[1]) & (Position <= region.coords[2]))
}

# In order to process each hmm file only once, our function for multi-threading needs to output all 
# objects needed from each hmm file
process_2R_regions <- function(sample.name){
	hmm.table <- load.hmm.file(sample.name, paste(coverage.folder, '2R/HMM_output', sep = '/'))
	output.list <- list(ace1 = shrink.compact.hmm(region.coords$ace1, hmm.table),
	                    cyp6 = shrink.compact.hmm(region.coords$cyp6, hmm.table)
	                   )
	return(output.list)
}

process_3R_regions <- function(sample.name){
	hmm.table <- load.hmm.file(sample.name, paste(coverage.folder, '3R/HMM_output', sep = '/'))
	output.list <- list(cyp6mz = shrink.compact.hmm(region.coords$cyp6mz, hmm.table),
	                    gste = shrink.compact.hmm(region.coords$gste, hmm.table)
	                   )
	return(output.list)
}

process_X_regions <- function(sample.name){
	hmm.table <- load.hmm.file(sample.name, paste(coverage.folder, 'X/HMM_output', sep = '/'))
	return(shrink.compact.hmm(region.coords$cyp9k1, hmm.table))
}

# Function that will calculate the mode of the HMM coverage in a given region
get.gene.mode <- function(hmm.data, target.region, window.size = 300){
	target.region <- as.numeric(target.region)
	hmm.in.region <- subset(hmm.data, Position >= (target.region[1] - window.size) & Position < target.region[2])$CNV
	# We calculate the mode. Where there is more than one mode, we take the smallest one
	un <- unique(hmm.in.region)
	tab <- tabulate(match(hmm.in.region, un))
	hmm.mode <- min(un[tab == max(tab)])
	hmm.mode
}

# Get the Agap numbers of the genes of interest, grouped by the cluster to which they belong.
focal.genes <- list(ace1   = c(Ace1 = 'AGAP001356'),
                    cyp6   = c(Cyp6aa1 = 'AGAP002862',
                               Cyp6aa2 = 'AGAP013128',
                               Coeae6o = 'AGAP002863',
                               AGAP002864 = 'AGAP002864',
                               Cyp6p1 = 'AGAP002868',
                               Cyp6p2 = 'AGAP002869',
                               Cyp6p3 = 'AGAP002865',
                               Cyp6p4 = 'AGAP002867',
                               Cyp6p5 = 'AGAP002866',
                               Cyp6ad1 = 'AGAP002870'),
                    cyp6mz = c(Cyp6m2 = 'AGAP008212',
                               Cyp6m3 = 'AGAP008213',
                               Cyp6m4 = 'AGAP008214',
                               Cyp6z1 = 'AGAP008219',
                               Cyp6z2 = 'AGAP008218',
                               Cyp6z3 = 'AGAP008217'),
                    gste   = c(Gste1 = 'AGAP009195',
                               Gste2 = 'AGAP009194',
                               Gste3 = 'AGAP009197',
                               Gste4 = 'AGAP009193',
                               Gste5 = 'AGAP009192',
                               Gste6 = 'AGAP009191',
                               Gste7 = 'AGAP009196',
                               Gstu4 = 'AGAP009190'),
                    cyp9k1 = c(Cyp9k1 = 'AGAP000818')
)

# Get the genetic coordinates of those genes.
gene.coords <- lapply(focal.genes, function(x) {X = all.gene.coordinates[x, ]; rownames(X) = names(x); X})

# Define the coordinates of the genetic regions that we will use for each cluster
region.coords = list(ace1   = c(3425000, 3650000),
                     cyp6   = c(28460000, 28800000),
                     cyp6mz = c(6900000, 7000000),
                     gste   = c(28570000, 28620000),
                     cyp9k1 = c(15220000, 15255000))

plotting.ranges <- list()
compact.hmm <- list()
cov.cnv.samples <- list()
# Create an empty dataframe where each row is a sample
hmm.cnv.table <- data.frame(row.names = sample.names)
# We will include a table giving the sex of the sample, since this is relevant to the Cyp9k1 copy number
hmm.cnv.table$sex <- meta[rownames(hmm.cnv.table), 'sex_call']
# And a column showing whether the sample has high coverage variance, since these will be dubious
hmm.cnv.table$high.var <- rownames(hmm.cnv.table) %in% high.var.samples

# Load up the 2R data and get the tables for the Ace1 and Cyp6 regions
cat('Creating shrunk table for Ace1 and Cyp6 regions.\n')
temp.2R <- mclapply(sample.names, process_2R_regions, mc.cores = num.cores)
# Because the above line uses parallel processing, I am not 100% confident that it will correctly preserve 
# sample order, so in the next couple of lines we use lapply on the sample.names vector to make sure the order 
# is correct. 
compact.hmm$ace1 <- lapply(sample.names, function(s) temp.2R[[s]]$ace1)
compact.hmm$cyp6 <- lapply(sample.names, function(s) temp.2R[[s]]$cyp6)
rm(temp.2R) 
# Get the ranges within this region that we wish to use for plotting
plotting.ranges$ace1 <- region.coords$ace1
plotting.ranges$cyp6 <- c(region.coords$cyp6[1], 28570000)
# Record the CNVs that were discovered based on coverage in these regions
hmm.cnv.table[['Ace1']] <- unlist(lapply(compact.hmm$ace1, get.gene.mode, target.region = gene.coords$ace1['Ace1', 2:3])) - 2
cov.cnv.samples$ace1 <- rownames(subset(hmm.cnv.table, Ace1 > 0))
cyp6aap.genes <- c('Cyp6aa1', 'Cyp6aa2', 'Cyp6p1', 'Cyp6p2', 'Cyp6p3', 'Cyp6p4', 'Cyp6p5')
for (g in cyp6aap.genes)
	hmm.cnv.table[[g]] <- unlist(lapply(compact.hmm$cyp6, get.gene.mode, target.region = gene.coords$cyp6[g, 2:3])) - 2
hmm.cnv.table[["Max_Cyp6aap"]] <- apply(hmm.cnv.table[, cyp6aap.genes], 1, max)
cov.cnv.samples$cyp6 <- rownames(subset(hmm.cnv.table, Max_Cyp6aap > 0))

# Load up the 3R data and get the tables for the Cyp6m-z and Gste regions
cat('Creating shrunk table for Cyp6M2-Z1 and Gste regions.\n')
temp.3R <- mclapply(sample.names, process_3R_regions, mc.cores = num.cores)
compact.hmm$cyp6mz <- lapply(sample.names, function(s) temp.3R[[s]]$cyp6mz)
compact.hmm$gste <- lapply(sample.names, function(s) temp.3R[[s]]$gste)
rm(temp.3R) 
plotting.ranges$cyp6mz <- region.coords$cyp6mz
plotting.ranges$gste <- region.coords$gste
for (g in names(focal.genes$cyp6mz))
	hmm.cnv.table[[g]] <- unlist(lapply(compact.hmm$cyp6mz, get.gene.mode, target.region = gene.coords$cyp6mz[g, 2:3])) - 2
hmm.cnv.table[["Max_Cyp6mz"]] <- apply(hmm.cnv.table[, names(focal.genes$cyp6mz)], 1, max)
cov.cnv.samples$cyp6mz <- rownames(subset(hmm.cnv.table, Max_Cyp6mz > 0))
for (g in names(focal.genes$gste))
	hmm.cnv.table[[g]] <- unlist(lapply(compact.hmm$gste, get.gene.mode, target.region = gene.coords$gste[g, 2:3])) - 2
hmm.cnv.table[["Max_Gstue"]] <- apply(hmm.cnv.table[, names(focal.genes$gste)], 1, max)
cov.cnv.samples$gste <- rownames(subset(hmm.cnv.table, Max_Gstue > 0))

# Load up the X data and get the tables for Cyp9k1	
cat('Creating shrunk table for Cyp9k1 region.\n')
temp.X <- mclapply(sample.names, process_X_regions, mc.cores = num.cores)
compact.hmm$cyp9k1 <- temp.X[sample.names]
rm(temp.X) 
plotting.ranges$cyp9k1 <- region.coords$cyp9k1
hmm.cnv.table[['Cyp9k1']] <- unlist(lapply(compact.hmm$cyp9k1, get.gene.mode, target.region = gene.coords$cyp9k1['Cyp9k1', 2:3])) - expected.copy.number.on.X
cov.cnv.samples$cyp9k1 <- rownames(subset(hmm.cnv.table, Cyp9k1 > 0))

# Get a names vector of the types of discordant reads that will be used 
discordant.read.types <- c(FA = 'FA', SS = 'SS', FM = 'FM', XC = 'crosschrom')

# Set the diagnostic reads that will be used to detect each CNV
known.cnvs <- list(ace1   = list(Dup1  = list(FA = matrix(c(3436850,3639550,3437150,3639850), 2, 2),
					                         #BP = data.frame(pos = c(3436926, 3639836), seq = c('GCGAA', 'GGAAT')))
                                              BP = data.frame(pos = c(3436926, 3639836), seq = c('GCGAA', 'TTGTT'))),
                                 Dup2  = list(FA = matrix(c(3447950,3518600,3448250,3518900), 2, 2)),
                                 Del1  = list(FM = matrix(c(3501850,3598750,3502150,3599050), 2, 2)),
                                 Del2  = list(FM = matrix(c(3539300,3573450,3539600,3573750), 2, 2)),
                                 Del3  = list(FM = matrix(c(3535850,3618700,3536150,3619000), 2, 2)),
                                 Del4  = list(FM = matrix(c(3512200,3615990,3512500,3616290), 2, 2))),
                   cyp6   = list(Dup1  = list(FA = matrix(c(28480150, 28483200, 28480450, 28483550), 2, 2)),
                                 Dup1a = list(BP = data.frame(pos = c(28480189, 28483475), 
                                                              seq = c('CGTAG', 'AATTG'))),
                                 Dup1b = list(BP = data.frame(pos = c(28480193, 28483675), 
                                                              seq = c('CTGCT', 'CCTTC'))),
                                 Dup2  = list(FA = matrix(c(28493450, 28497000, 28493750, 28497300), 2, 2),
                                              BP = data.frame(pos = c(28493547, 28497279), 
                                                              seq = c('GCCGC','TTTAA'))),
                                 Dup3  = list(FA = matrix(c(28479350, 28483100, 28479650, 28483400), 2, 2),
                                              BP = data.frame(pos = c(28479407, 28483372),
                                                              seq = c('GCTTA', 'CAAAG'))),
                                 Dup4  = list(FA = matrix(c(28478850, 28482750, 28479150, 28483050), 2, 2),
                                              BP = data.frame(pos = c(28478925, 28483069),
                                                              seq = c('TACTT', 'CATGT'))),
                                 Dup5  = list(FA = matrix(c(28480300, 28484200, 28480600, 28484500), 2, 2), 
                                              BP = data.frame(pos = c(28480372, 28484518),
                                                              seq = c('AAGAG', 'ACAAA'))),
                                 Dup6  = list(FA = matrix(c(28478150, 28483850, 28478450, 28484150), 2, 2),
                                              BP = data.frame(pos = c(28478272, 28484157),
                                                              seq = c('ATCAC', 'CTAGA'))),
                                 Dup7  = list(SS = matrix(c(28478000, 28486000, 28478300, 28486300), 2, 2),
                                              BP = data.frame(pos = c(28478057, 28486036),
                                                              seq = c('AGAGC','TTTTT'))),
                                 Dup8  = list(FA = matrix(c(28475900, 28484700, 28476200, 28485000), 2, 2),
                                              BP = data.frame(pos = c(28475996, 28485005),
                                                              seq = c('AGCGA', 'CAAAT'))),
                                 Dup9  = list(FA = matrix(c(28479100, 28491200, 28479400, 28491500), 2, 2),
                                              BP = data.frame(pos = c(28479181, 28491431),
                                                              seq = c('TGTTC', 'TGTGG'))),
                                 Dup10 = list(FA = matrix(c(28477800, 28490850, 28478100, 28491150), 2, 2),
                                              BP = data.frame(pos = c(28477889, 28491215),
         # It turns out that Dup10 is not quite as simple as made out in the Supp Mat for the Genome Research paper.
         # There is actually some kind of insertion / mutation happening around the breakpoint, and the aligners used
         # for this experiment and in Ag1000G deal with this differently (so we need to check the phase3 data and 
         # onwards to see what happens there since they used a different aligner to phase 2). We therefore need to change
         # the sequence slightly here.
                                                              seq = c('TGTAG','AACTT'))),
                                                             #seq = c('TGTAG','ACTCT'))),
                                 Dup11 = list(FA = matrix(c(28487450, 28517800, 28487750, 28518100), 2, 2),
                                              BP = data.frame(pos = c(28487546, 28518123),
                                                              seq = c('AACAC', 'TTATC'))),
                                 Dup12 = list(FA = matrix(c(28474450, 28519650, 28474750, 28519950), 2, 2),
                                              BP = data.frame(pos = c(28474576, 28520016),
                                                              seq = c('CCGAC', 'ACGGT'))),
                                 Dup13 = list(FA = matrix(c(28472650, 28522350, 28472950, 28522650), 2, 2),
                                              BP = data.frame(pos = c(28472728, 28522671),
                                                              seq = c('ACCGC', 'AGCTG'))),
                                 Dup14 = list(FA = matrix(c(28473800, 28563200, 28474100, 28563500), 2, 2),
                                              BP = data.frame(pos = c(28473874, 28563596),
                                                              seq = c('CCCAC', 'AGTTG'))),
                                 Dup15 = list(FA = matrix(c(28465600, 55958800, 28465900, 55959100), 2, 2),
                                              BP = data.frame(pos = c(28465673, NA),
                                                              seq = c('CAGCC', NA))),
                                 Dup16 = list(FA = matrix(c(28480500, 28483300, 28480800, 28483600), 2, 2),
                                              BP = data.frame(pos = c(28480547, 28483236),
                                                               seq = c('CCATT', 'TTAGT'))),
                                 Dup17 = list(FA = matrix(c(28477500, 28484900, 28477800, 28485200), 2, 2),
                                              BP = data.frame(pos = c(28477540, 28485380),
                                                              seq = c('TGCTG', 'ATCGG'))),
                                 Dup18 = list(FA = matrix(c(28479500, 28494200, 28479800, 28494500), 2, 2),
                                              BP = data.frame(pos = c(28479548, 28494597),
                                                              seq = c('AGTCG', 'TTGTC'))),
                                 Dup19 = list(FA = matrix(c(28475480, 28556250, 28475780, 28556550), 2, 2),
                                              BP = data.frame(pos = c(28475490, 28556726),
                                                              seq = c('AATAG', 'TGTGT'))),
                                 Dup20 = list(FA = matrix(c(28473590, 28794750, 28473890, 28795050), 2, 2),
                                              BP = data.frame(pos = c(28473600, 28795255),
                                                              seq = c('ATACT', 'CAAAA'))),
                                 Dup21 = list(FA = matrix(c(28475100, 28483000, 28475400, 28483300), 2, 2),
                                              BP = data.frame(pos = c(28475128, 28483473),
                                                              seq = c('AGCCG', 'TGCAC'))),
                                 Dup22 = list(FA = matrix(c(28477200, 28484000, 28477500, 28484300), 2, 2),
                                              BP = data.frame(pos = c(28477220, 28484338),
                                                              seq = c('GTGGA', 'CGAGT'))),
                                 Dup23 = list(FA = matrix(c(28497300, 28371800, 28497600, 28372100), 2, 2),
                                              BP = data.frame(pos = c(NA, 28497740),
                                                              seq = c(NA, 'TTGGC'))),
                                 Dup24  = list(SS = matrix(c(28479500, 28483000, 28479800, 28483300), 2, 2),
                                               BP = data.frame(pos = c(28480585, 28483442),
                                                               seq = c('AAACA','TTAAC'))),
				   	             # For 25 and 26, the FA reads would overlap, so we just use the BP reads. 
                                 Dup25 = list(BP = data.frame(pos = c(28480335, 28483384),
                                                              seq = c('GGCGT', 'CATAT'))),
                                 Dup26 = list(BP = data.frame(pos = c(28480166, 28483253),
                                                              seq = c('AACGT', 'TGTGT'))),
                                 Dup27 = list(FA = matrix(c(28496700, 28498900, 28497000, 28499200), 2, 2)),
                                 Dup28 = list(FA = matrix(c(28477700, 28496600, 28478000, 28496900), 2, 2),
                                              BP = data.frame(pos = c(28477710, 28496953),
                                                              seq = c('CTGTA', 'ATTCT'))),
                                 Dup29 = list(FA = matrix(c(28494000, 28496000, 28494300, 28496300), 2, 2),
                                              BP = data.frame(pos = c(28494017, 28496505),
                                                              seq = c('TGGAA', 'TTTGC'))),
                                 Dup30 = list(FA = matrix(c(28478900, 28484700, 28479200, 28485000), 2, 2),
                                              BP = data.frame(pos = c(28478987, 28485033),
                                                              seq = c('AACAG', 'ACGTT')))),
                   cyp6mz = list(Dupm1 = list(BP = data.frame(pos = c(6927942, NA),
                                                              seq = c('ATTAT', NA))),
                                 Dupm2 = list(FA = matrix(c(6933100, 6934900, 6933400, 6935200), 2, 2)),
                                 Dupm3 = list(FA = matrix(c(6929600, 6932500, 6929900, 6932800), 2, 2)),
                                 Dupm4 = list(FA = matrix(c(6929900, 6936600, 6930200, 6936900), 2, 2),
                                              BP = data.frame(pos = c(6929933, 6936902),
                                                              seq = c('TTAAA', 'TGTCG'))),
                                 Dupm5 = list(FA = matrix(c(6933800, 6938300, 6934100, 6938600), 2, 2),
                                              BP = data.frame(pos = c(6933972, 6938659),
                                                              seq = c('AAACC', 'GTCGG'))),
                                 Dupz1 = list(FA = matrix(c(6968950, 6979300, 6969250, 6979600), 2, 2),
                                              BP = data.frame(pos = c(6968962, 6979681),
                                                              seq = c('ACGCT', 'AGGTT'))),
                                 Dupz2 = list(FA = matrix(c(6975100, 6977100, 6975400, 6977400), 2, 2),
                                              BP = data.frame(pos = c(NA, 6977514), # clipped reads align at 6975066
                                                              seq = c(NA, 'TAAGA'))),
                                 Dupz3 = list(FA = matrix(c(6971450, 6977800, 6971750, 6978100), 2, 2),
                                              BP = data.frame(pos = c(6971484, NA),  # clipped reads align at 6978084
                                                              seq = c('GCAAA', NA))),
                                 Dupz4 = list(FA = matrix(c(6972700, 6977350, 6973000, 6977650), 2, 2),
                                              BP = data.frame(pos = c(6972775, 6977699),  
                                                              seq = c('GAATG', 'GTCCA'))),
                                 Dupz5 = list(FA = matrix(c(6969800, 6975700, 6970100, 6976000), 2, 2)),
                                 Dupmz1 = list(FA = matrix(c(6982700, 6879900, 6983000, 6880200), 2, 2))),
                   gste   = list(Dup1  = list(FA = matrix(c(6968950, 6979300, 6969250, 6979600), 2, 2),
                                              BP = data.frame(pos = c(28596818, 28598850),
                                                              seq = c('TTTTG', 'CGTTT'))),
                                 # The following definition includes some false positives. 
                                 Dup2  = list(BP.weak = data.frame(pos = c(28596390, 28598923),
                                                                   seq = c('GGGGG', 'TTCCC'))),
                                 Dup3  = list(FA = matrix(c(28590500, 28592950, 28590800, 28593250), 2, 2),
                                              BP = data.frame(pos = c(28590597, 28593254),
                                                              seq = c('TCAAA', 'AGGGC'))),
                                 Dup4  = list(FA = matrix(c(28595050, 28598750, 28595350, 28599050), 2, 2),
                                              BP = data.frame(pos = c(28595162, 28599081),
                                                              seq = c('TTCTA', 'AGAAC'))),
                                 Dup5  = list(BP = data.frame(pos = c(28593122, 28598971),
                                                              seq = c('GTCAT', 'ATTTA'))),
                                 Dup6  = list(FA = matrix(c(28596250, 28601900, 28596550, 28602200), 2, 2),
                                              BP = data.frame(pos = c(28596241, 28602177),
                                                              seq = c('ACAAC', 'GAAGC'))),
                                 # For the XC reads, the first two rows are the CNV start point, and the second two rows are the end point
                                 Dup7  = list(XC = data.frame(c(28597400, 3696450, 28603950, 26597300), c(28597700, 3696750, 28604250, 26597600), c('3R', '2L', '3R', 'UNKN')),
                                              BP = data.frame(pos = c(28597504, 28604250),
                                                              seq = c('GTCCA', 'GCTGT'))),
                                 Dup8  = list(BP = data.frame(pos = c(28594797, 28602349),
                                                              seq = c('GTCCC', 'CAGGG'))),
                                 Dup9  = list(FA = matrix(c(28591050, 28600850, 28591350, 28601150), 2, 2),
                                              BP = data.frame(pos = c(28591140, 28601188),
                                                              seq = c('AGAAG', 'GATGA'))),
                                 Dup10 = list(FA = matrix(c(28593550, 28603350, 28593850, 28603650), 2, 2),
                                              BP = data.frame(pos = c(28593642, 28603786),
                                                              seq = c('TCGCT', 'AAGAC'))),
                                 Dup11 = list(XC = data.frame(c(28581250, 29210650, 28604650, 29210650), c(28581550, 29210950, 28604950, 29210950), c('3R', 'UNKN', '3R', 'UNKN')),
                                              BP = data.frame(pos = c(28581256, 28604994),
                                                              seq = c('CCATT', 'GGTAA'))),
                                 Dup12 = list(FA = matrix(c(28597000, 28599950, 28597300, 28600250), 2, 2),
                                              BP = data.frame(pos = c(28597030, 28600292),
                                                              seq = c('TACTG', 'CATCT'))),
                                 Dup13 = list(FA = matrix(c(28597000, 28598900, 28597300, 28599200), 2, 2),
                                              BP = data.frame(pos = c(28597181, NA), # clipped sequences align at 28599287
                                                              seq = c('TACTC', NA))),
                                 Dup14 = list(FA = matrix(c(28599800, 28607200, 28600100, 28607500), 2, 2),
                                              BP = data.frame(pos = c(28599926, 28607500), 
                                                              seq = c('CGACG', 'ATGCA'))),
                                 Dup15 = list(FA = matrix(c(28596200, 28598500, 28596500, 28598800), 2, 2),
                                              BP = data.frame(pos = c(28596252, 28598948),
                                                              seq = c('TTGGA', 'TTGAC'))),
                                 Dup16 = list(FA = matrix(c(28597300, 28603200, 28597600, 28603500), 2, 2),
                                              BP = data.frame(pos = c(28597383, 28603517), 
                                                              seq = c('ACATT', 'ATTAC')))),
                   cyp9k1 = list(Dup1  = list(FA = matrix(c(15242500, 15244500, 15242800, 15244800), 2, 2),
                                              BP = data.frame(pos = c(15242505,15244812),
                                                            seq = c('GTTTG', 'CATAT'))),
                                 Dup2  = list(FA = matrix(c(15238300, 15240800, 15238600, 15241100), 2, 2),
                                              BP = data.frame(pos = c(15238400, 15241082),
                                                              seq = c('CCGGC',' CGGTA'))),
                                 Dup3  = list(FA = matrix(c(15240300, 15243450, 15240600, 15243750), 2, 2),
                                              BP = data.frame(pos = c(NA, 15243860),
                                                              seq = c(NA, 'TGAAC'))),
                                 Dup4  = list(FA = matrix(c(15240600, 15244200, 15240900, 15244500), 2, 2),
                                              BP = data.frame(pos = c(15240608, 15244503),
                                                              seq = c('ATAAA', 'ACTGG'))),
                                 Dup5  = list(FA = matrix(c(15238800, 15243850, 15239100, 15244150), 2, 2),
                                              BP = data.frame(pos = c(15238911, 15244175),
                                                              seq = c('CACGT', 'AGTAA'))),
                                 Dup6  = list(FA = matrix(c(15236400, 15243250, 15236700, 15243550), 2, 2),
                                              BP = data.frame(pos = c(15236449, 15243646),
                                                              seq = c('TTTTT', 'GTTTT'))),
                                 Dup7  = list(SS = matrix(c(15245400, 15246900, 15245700, 15247200), 2, 2),
                                              BP = data.frame(pos = c(15245768, 15247258),
                                                              seq = c('TTTGT', 'TCTAA'))),
                                 Dup8  = list(FA = matrix(c(15239200, 15247250, 15239500, 15247550), 2, 2),
                                              BP = data.frame(pos = c(15239276, 15247645),
                                                              seq = c('AACAT', 'TTGCT'))),
                                 Dup9  = list(FA = matrix(c(15239100, 15248900, 15239400, 15249200), 2, 2),
                                              BP = data.frame(pos = c(15239184, 15249314),
                                                              seq = c('GCACA', 'AGTAC'))),
                                 Dup10 = list(FA = matrix(c(15234900, 15244750, 15235200, 15245050), 2, 2),
                                              BP = data.frame(pos = c(15234989, 15245128),
                                                              seq = c('GCACC', 'ATTCT'))),
                                 Dup11 = list(FA = matrix(c(15236900, 15246800, 15237200, 15247100), 2, 2),
                                              BP = data.frame(pos = c(15236922, 15247159),
                                                              seq = c('CATTA', 'TATCT'))),
                                 Dup12 = list(FA = matrix(c(15234400, 15244350, 15234700, 15244650), 2, 2),
                                              BP = data.frame(pos = c(15234434, 15244702),
                                                              seq = c('AACAG', 'TACTA'))),
                                 Dup13 = list(FA = matrix(c(15240100, 15250250, 15240400, 15250550), 2, 2),
                                              BP = data.frame(pos = c(15240067, 15250575),
                                                              seq = c('CCTAA', 'GTGTA'))),
                                 # Dup 14 seems to have a different endpoint to Dup15, but the same insertion point way upstream
                                 Dup14 = list(FA = matrix(c(15244200, 9676400, 15244500, 9676700), 2, 2),
                                             # Because Dup14 has the same start pos as Dup15, we don't use the start breakpoint as a diagnostic
                                             #BP = data.frame(pos = c(15233807, 15244936),
                                             #                seq = c('GGGTT', 'CCCAA'))),
                                              BP = data.frame(pos = c(NA, 15244936),
                                                              seq = c(NA, 'CCCAA'))),
                                 Dup15 = list(FA = matrix(c(15246250, 9676400, 15246550, 9676700), 2, 2),
                                             # Because Dup14 has the same start pos as Dup15, we don't use the start breakpoint as a diagnostic
                                             #BP = data.frame(pos = c(15233807, 15246640),
                                             #                seq = c('GGGTT', 'CCCAA'))),
                                              BP = data.frame(pos = c(NA, 15246640),
                                                              seq = c(NA, 'CCCAA'))),
                                 Dup16 = list(FA = matrix(c(15222700, 15244300, 15223000, 15244600), 2, 2),
                                              BP = data.frame(pos = c(NA, 15244755),
                                                              seq = c(NA, 'AAGTA'))),
                                 Dup17 = list(FA = matrix(c(15237150, 15243650, 15237450, 15243950), 2, 2),
                                              BP = data.frame(pos = c(15237138, 15243975),
                                                              seq = c('TTGCT', 'TTTCG'))),
                                 Dup18 = list(FA = matrix(c(15236100, 15243500, 15236400, 15243800), 2, 2),
                                              BP = data.frame(pos = c(NA, 15243915),  # clipped sequence aligns to 15236175
                                                              seq = c(NA, 'CGGCG'))),
                                 Dup19 = list(FA = matrix(c(15238800, 15251100, 15239100, 15251400), 2, 2),
                                              BP = data.frame(pos = c(15238878, 15251503),  
                                                              seq = c('TAAAT', 'GTTAC'))),
                                 Dup20 = list(FA = matrix(c(15237350, 15243100, 15237650, 15243400), 2, 2),
                                              BP = data.frame(pos = c(15237397, 15243514),  
                                                              seq = c('ATGTT', 'TTACG'))),
                                 Dup21 = list(FA = matrix(c(15237450, 15250300, 15237750, 15250600), 2, 2),
                                              BP = data.frame(pos = c(15237482, 15250699),  
                                                              seq = c('CTCTG', 'TTCTC'))),
                                 Dup22 = list(FA = matrix(c(15240650, 15250300, 15240950, 15250600), 2, 2),
                                              BP = data.frame(pos = c(15240680, 15250670),  
                                                              seq = c('TTCCA', 'ATTCT'))),
                                 Dup23 = list(FA = matrix(c(15241800, 15248100, 15242100, 15248400), 2, 2),
                                              BP = data.frame(pos = c(15241929, 15248352),  
                                                              seq = c('AACAA', 'CACGT'))),
                                 Dup24 = list(FA = matrix(c(15238550, 15254800, 15238850, 15255100), 2, 2)),
                                 Dup25 = list(FA = matrix(c(15223500, 15246350, 15223800, 15246650), 2, 2)),
                                 Dup26 = list(FA = matrix(c(15222700, 15247750, 15223000, 15248050), 2, 2)),
                                 Dup27 = list(FA = matrix(c(15237400, 15248350, 15237700, 15248650), 2, 2),
                                              BP = data.frame(pos = c(15237566, NA),  
                                                              seq = c('AATGT', NA))),
                                 Dup28 = list(BP = data.frame(pos = c(NA, 15246640),  
                                                              seq = c(NA, 'TCGAG'))))
)

# Write a function to load the results of discordant read analysis from file for a given sample
get.discordant.reads <- function(this.sample, target.folder, target.region.coords){
	this.file <- grep(paste(this.sample, '.*csv$', sep = ''), list.files(target.folder), value = T)
	if (length(this.file) != 1)
		stop('There should be one and only one file in the folder that matches the sample name.')
   	pos.table <- read.table(paste(target.folder, this.file, sep='/'), header = T, sep='\t', colClasses = c('character', 'integer', 'character'))
	split.pos.table <- split(pos.table, pos.table$Type)
	get.discordant.read.type <- function(read.type){
		full.read.type <- paste(read.type, 'mapq >= 10')
		if (!(full.read.type %in% names(split.pos.table))){
			return(data.frame(Position = integer(), Mate.position = integer()))
		}
		allpos <- split.pos.table[[full.read.type]][, c('Position', 'Mate.position')]
		if (read.type %in% c('FA', 'SS', 'FM')){
			allpos$Mate.position <- as.integer(allpos$Mate.position)
			which.in.region <- ((allpos$Position >= target.region.coords[1]) & (allpos$Position <= target.region.coords[2])) | ((allpos$Mate.position >= target.region.coords[1]) & (allpos$Mate.position <= target.region.coords[2]))
		}
		else{
			which.in.region <- (allpos$Position >= target.region.coords[1]) & (allpos$Position <= target.region.coords[2])
		}
		allpos <- allpos[which.in.region, ]
		allpos[order(allpos$Position),]
	}
	lapply(discordant.read.types, get.discordant.read.type)
}

# Write a function to load the results of breakpoint read analysis from file for a given sample
get.breakpoint.reads <- function(this.sample, target.folder, target.region.coords){
	this.file <- grep(paste(this.sample, '.*csv$', sep = ''), list.files(target.folder), value = T)
	if (length(this.file) != 1)
		stop('There should be one and only one file in the folder that matches the sample name.')
   	pos.table <- read.table(paste(target.folder, this.file, sep='/'), header = T, sep='\t', colClasses = c('character', 'integer', 'character'))
	which.in.region <- (pos.table$Position >= target.region.coords[1]) & (pos.table$Position <= target.region.coords[2])
	pos.table <- pos.table[which.in.region, ]
	clipping.start.point <- pos.table[pos.table$Type == 'soft clipping start point mapq >= 10', c('Position', 'Clipped_sequence')]
	clipping.end.point <- pos.table[pos.table$Type == 'soft clipping end point mapq >= 10', c('Position', 'Clipped_sequence')]
	#
	list(CSP = clipping.start.point[order(clipping.start.point$Position), ],
	     CEP = clipping.end.point[order(clipping.end.point$Position), ])
}

get.diagnostic.reads <- function(this.sample, disc.target.folder, bp.target.folder, target.region.coords){
	cat('\tLoading diagnostic reads for ', this.sample, '.\n', sep = '')
	output <- get.discordant.reads(this.sample, disc.target.folder, target.region.coords)
	output[['BP']] <- get.breakpoint.reads (this.sample, bp.target.folder, target.region.coords)
	output
}

# Write a function to load the results of discordant and breakpoint reads from file for a vector of samples
get.diagnostic.reads.allsamples <- function(sample.names, disc.target.folder, bp.target.folder, target.region.coords){
	cat('Loading discordant reads from ', disc.target.folder, ' and breakpoint reads from ', bp.target.folder, '.\n', sep = '')
	lapply(sample.names, get.diagnostic.reads, disc.target.folder, bp.target.folder, target.region.coords)
}

# Function to count the diagnostic reads for a given CNV. 
# di is a table of observed discordant and breakpoint reads for a given sample. 
count.reads.per.dup <- function(Dup.diagnostics, di){
	# There shouldn't be more than one of FA, SS, FM and XC
	if (sum(names(Dup.diagnostics) %in% c('FA', 'SS', 'FM', 'XC')) > 1)
		stop('There should not be more than one type of diagnostic read for a given Duplication.')
	num.reads <- setNames(numeric(7), c('FA', 'SS', 'FM', 'XCS', 'XCE', 'CSP', 'CEP'))
	if ('FA' %in% names(Dup.diagnostics)){
		num.reads['FA'] <- sum((di$FA$Position > Dup.diagnostics$FA[1,1]) & di$FA$Position < Dup.diagnostics$FA[1,2]
							 & (di$FA$Mate.position > Dup.diagnostics$FA[2,1]) & (di$FA$Mate.position < Dup.diagnostics$FA[2,2]))
	}
	if ('SS' %in% names(Dup.diagnostics)){
		num.reads['SS'] <- sum((di$SS$Position > Dup.diagnostics$SS[1,1]) & di$SS$Position < Dup.diagnostics$SS[1,2]
		                     & (di$SS$Mate.position > Dup.diagnostics$SS[2,1]) & (di$SS$Mate.position < Dup.diagnostics$SS[2,2]))
	}
	if ('FM' %in% names(Dup.diagnostics)){
		num.reads['FM'] <- sum((di$FM$Position > Dup.diagnostics$FM[1,1]) & di$FM$Position < Dup.diagnostics$FM[1,2]
							 & (di$FM$Mate.position > Dup.diagnostics$FM[2,1]) & (di$FM$Mate.position < Dup.diagnostics$FM[2,2]))
	}
	if ('XC' %in% names(Dup.diagnostics)){
		these.XC <- di$XC
		these.XC$Matechrom <- sub(':.*', '', these.XC$Mate.position)
		these.XC$Matepos <- as.integer(sub('.*:' , '', these.XC$Mate.position))
		num.reads['XCS'] <- sum((these.XC$Position > Dup.diagnostics$XC[1,1]) & (these.XC$Position < Dup.diagnostics$XC[1,2])
		                      & (these.XC$Matechrom == Dup.diagnostics$XC[2,3]) & (these.XC$Matepos > Dup.diagnostics$XC[2,1]) & (these.XC$Matepos < Dup.diagnostics$XC[2,2]), na.rm = T)
		num.reads['XCE'] <- sum((these.XC$Position > Dup.diagnostics$XC[3,1]) & (these.XC$Position < Dup.diagnostics$XC[3,2])
		                      & (these.XC$Matechrom == Dup.diagnostics$XC[4,3]) & (these.XC$Matepos > Dup.diagnostics$XC[4,1]) & (these.XC$Matepos < Dup.diagnostics$XC[4,2]), na.rm = T)
	}
	if ('BP' %in% names(Dup.diagnostics)){
		CEP.seq <- subset(di$BP$CEP, Position == Dup.diagnostics$BP$pos[1])$Clipped_sequence
		num.reads['CEP'] <- sum(substr(reverse(CEP.seq), 1, 5) == Dup.diagnostics$BP$seq[1])
		CSP.seq <- subset(di$BP$CSP, Position == Dup.diagnostics$BP$pos[2])$Clipped_sequence
		num.reads['CSP'] <- sum(substr(CSP.seq, 1, 5) == Dup.diagnostics$BP$seq[2])
	}
	# There is a special case in Gste2 where one Dup isn't defined very effectively by diagnostic
	# reads. The best diagnostic based on phase 2 data still includes some false positives, which 
	# were minimised by taking the minimum of CEP and CSP, rather than the sum
	if ('BP.weak' %in% names(Dup.diagnostics)){
		CEP.seq <- subset(di$BP$CEP, Position == Dup.diagnostics$BP$pos[1])$Clipped_sequence
		num.reads.cep <- sum(substr(reverse(CEP.seq), 1, 5) == Dup.diagnostics$BP$seq[1])
		CSP.seq <- subset(di$BP$CSP, Position == Dup.diagnostics$BP$pos[2])$Clipped_sequence
		num.reads.csp <- sum(substr(CSP.seq, 1, 5) == Dup.diagnostics$BP$seq[2])
		num.reads['BP.weak'] <- min(num.reads.cep, num.reads.csp)
	}
	num.reads
}

# Function to apply count.reads.per.dup to all CNVs for a given sample
count.diagnostic.reads <- function(diagnostic.reads.list, known.cnvs.list, all.breakpoint.positions, verbose = F){
	# Subset the BP reads so that they only include reads that match at least one of
	# the breakpoints. That will greatly reduce the list that needs searching for every
	# CNV. 
	diagnostic.reads.list$BP$CSP <- subset(diagnostic.reads.list$BP$CSP, Position %in% all.breakpoint.positions)
	diagnostic.reads.list$BP$CEP <- subset(diagnostic.reads.list$BP$CEP, Position %in% all.breakpoint.positions)
	num.reads <- lapply(known.cnvs.list, count.reads.per.dup, diagnostic.reads.list)
	if (verbose)
		print(num.reads)
	sapply(num.reads, sum)
}

# Function to apply counts.diagnostic.reads to all samples
count.diagnostic.reads.allsamples <- function(diagnostic.reads.allsamples, known.cnvs.list){
	# Create a single object containing all of the breakpoint positions
	breakpoint.summary <- unlist(sapply(known.cnvs.list, function(x) x$BP$pos))
	t(sapply(diagnostic.reads.allsamples, count.diagnostic.reads, known.cnvs.list, breakpoint.summary))
}

# Function to produce a table of all observed CNVs in all samples, given a list of coverage-based calls,
# a table of counts of diagnostic reads for each CNV in all samples, and a threshold number of diagnostic
# reads required to call a CNV
get.read.based.cnvs <- function(coverage.cnv.calls, diagnostic.read.counts.table, threshold.diagnostic.reads = 4){
	read.based <- cbind(diagnostic.read.counts.table >= threshold.diagnostic.reads)
	read.cnv.samples <- rownames(read.based)[apply(read.based, 1, any)]
	# Unlike in phase2, we use Dup0 to indicate that there is ANY coverage call, not just the samples that
	# have a coverage call but not a read-based call
	cbind(Dup0 = rownames(read.based) %in% coverage.cnv.calls, read.based)
}

# Set the folders in which to look for FA, SS, FM and XC reads
SSFA.folders <- setNames(paste(diagnostic.reads.folder, c('SSFA/2R/Ace1_region',
                                                          'SSFA/2R/Cyp6_region',
                                                          'SSFA/3R/Cyp6zm_region',
                                                          'SSFA/3R/Gste_region',
                                                          'SSFA/X/Cyp9k1_region'),
                                                        sep = '/'),
                         c('ace1', 'cyp6', 'cyp6mz', 'gste', 'cyp9k1')
                        )

# Set the folders in which to look for breakpoint reads. 
breakpoints.folders <- setNames(paste(diagnostic.reads.folder, c('breakpoints/2R/Ace1_region',
                                                                 'breakpoints/2R/Cyp6_region',
                                                                 'breakpoints/3R/Cyp6zm_region',
                                                                 'breakpoints/3R/Gste_region',
                                                                 'breakpoints/X/Cyp9k1_region'),
                                                               sep = '/'),
                                c('ace1', 'cyp6', 'cyp6mz', 'gste', 'cyp9k1')
                               )

cat('Loading discordant and breakpoint reads\n')
diagnostic.reads <- mcmapply(get.diagnostic.reads.allsamples, SSFA.folders, breakpoints.folders, region.coords, MoreArgs = list(sample.names = sample.names), SIMPLIFY = F, mc.preschedule = F, mc.cores = num.cores)

cat('Counting diagnostic reads\n')
diagnostic.read.counts <- mcmapply(count.diagnostic.reads.allsamples, diagnostic.reads, known.cnvs, mc.preschedule = F, mc.cores = num.cores)
cat('Calling read-based CNVs\n')
read.based.cnvs <- mapply(get.read.based.cnvs, cov.cnv.samples, diagnostic.read.counts, SIMPLIFY = F)
# As a special case for cyp6, we don't score the presence of Dup1 if either Dup1a or Dup1b are present (because Dup1
# just represents an ambiguous call for Dup1a / Dup1b. 
cyp6.Dup1ab <- read.based.cnvs$cyp6[, 'Dup1a'] | read.based.cnvs$cyp6[, 'Dup1b']
read.based.cnvs$cyp6[cyp6.Dup1ab, 'Dup1'] <- F

# Combine the separate CNV tables into a single table
full.cnv.table <- do.call(cbind, read.based.cnvs)
full.cnv.table <- cbind(rownames(full.cnv.table) %in% high.var.samples, full.cnv.table)
gene.cluster.names <- setNames(c('Ace1', 'Cyp6aap', 'Cyp6mz', 'Gstue', 'Cyp9k1'), 
                               names(diagnostic.read.counts))
# The column names for the full table will be with the gene cluster codes that we used in the phase 2 analysis,
# and will split the Cyp6mz cnvs into their constituent clusters
colnames(full.cnv.table) <- c('High.var.sample', stri_replace_all_fixed(unlist(sapply(names(gene.cluster.names), function(x) paste(gene.cluster.names[x], colnames(read.based.cnvs[[x]]), sep = '_'))),
                                                                        c('Cyp6mz_Dupm[^z]', 'Cyp6mz_Dupz'),
                                                                        c('Cyp6m_Dup', 'Cyp6z_Dup'),
                                                                        vectorize_all = F)
                             )
# Get the results in a different format. Here as a list where, for each CNV, we have a vector of samples that 
# carry it. 
cnv.sample.lists <- lapply(read.based.cnvs, function(m) apply(m, 2, function(x) rownames(m)[x]))

# Get a vector of samples that carry at least one CNV based on read calls.
# Dup0 is always the first column, so the -1 index removes Dup0
read.cnv.samples <- lapply(read.based.cnvs, function(x) rownames(x)[apply(x[, -1], 1, any)])

# For each cluster, get a vector of samples that have been called as carrying a cnv by reads but not coverage,
# ignoring samples with high variance.
cov.based.negatives <- lapply(mapply(setdiff, read.cnv.samples, cov.cnv.samples), setdiff, high.var.samples)
# And vice-versa
read.based.negatives <- lapply(mapply(setdiff, cov.cnv.samples, read.cnv.samples), setdiff, high.var.samples)

dir.create('target_regions_analysis', showWarnings = FALSE)
write.table(full.cnv.table, file = 'target_regions_analysis/focal_region_CNV_table.csv', sep = '\t', col.names = NA, quote = F)
write.table(hmm.cnv.table, file = 'target_regions_analysis/HMM_gene_copy_number.csv', sep = '\t', col.names = NA, quote = F)

source(plotting.functions.file)
save.image('target_regions_analysis/target_regions_analysis.Rdata')

