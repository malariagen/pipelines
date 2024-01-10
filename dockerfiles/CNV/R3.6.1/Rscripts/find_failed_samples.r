arg.values <- commandArgs(trailingOnly=T)
cat('Arguments:', arg.values, '\n', sep = '\n')

study.id <- arg.values[1]
bampaths.table.filename <- arg.values[2]
expected.coverage.files <- as.numeric(arg.values[3])
expected.diagnostic.reads.files <- as.numeric(arg.values[4])
output.filename <- arg.values[5]

working.folder <- paste('~/personal', study.id, sep = '/')

bampaths.table <- read.table(bampaths.table.filename, sep = '\t', header = F, row.names = 1)
expected.samples <- setNames(nm = rownames(bampaths.table))

coverage.folder <- paste(working.folder, 'coverage', sep = '/')
all.coverage.files <- list.files(coverage.folder, recursive = T)
csv.coverage.files <- all.coverage.files[grepl('.csv$', all.coverage.files)]
matching.coverage.files <- sapply(expected.samples, function(x) sum(grepl(x, csv.coverage.files)))
missing.coverage.samples <- matching.coverage.files[matching.coverage.files < expected.coverage.files]

cat('\nThe following samples did not have the full complement of coverage output files (', 
     expected.coverage.files, ' files were expected per sample, the number observed for ',
    'each failing sample is shown below).\n\n', sep = '')
print(missing.coverage.samples)

diagnostic.reads.folder <- paste(working.folder, 'diagnostic_reads', sep = '/')
all.diagnostic.reads.files <- list.files(diagnostic.reads.folder, recursive = T)
csv.diagnostic.reads.files <- all.diagnostic.reads.files[grepl('.csv$', all.diagnostic.reads.files)]
matching.diagnostic.reads.files <- sapply(expected.samples, function(x) sum(grepl(x, csv.diagnostic.reads.files)))
missing.diagnostic.reads.samples <- matching.diagnostic.reads.files[matching.diagnostic.reads.files < expected.diagnostic.reads.files]

cat('\n\nThe following samples did not have the full complement of diagnostic read output files (', 
     expected.diagnostic.reads.files, ' files were expected per sample, the number observed for ',
    'each failing sample is shown below).\n\n', sep = '')
print(missing.diagnostic.reads.samples)

all.missing.samples <- sort(unique(c(names(missing.coverage.samples), names(missing.diagnostic.reads.samples))))
cat('\n\nOverall,', length(all.missing.samples), 'samples were missing at least one output file.\n')

missing.samples.bampaths <- bampaths.table[all.missing.samples, , drop = F]
write.table(missing.samples.bampaths, output.filename, sep = '\t', col.names = F, row.names = T, quote = F)

