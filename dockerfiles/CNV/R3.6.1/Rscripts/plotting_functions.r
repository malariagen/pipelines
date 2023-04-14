library(RColorBrewer)

# Create a vector of colours to use for plotting. 
# The last two CNVs are named Dup1a and Dup1b, for the Cyp6 cluster
n.colours <- 30
duplication.colours <- c(colorRampPalette(brewer.pal(12, 'Paired'))(n.colours), 'grey30', 'grey70')
#duplication.colours <- c(rgb(0.6509804, 0.8078431, 0.8901961), rgb(0.1215686, 0.4705882, 0.7058824), rgb(0.6980392, 0.8745098, 0.5411765), rgb(0.2, 0.627451, 0.172549), rgb(0.9843137, 0.6039216, 0.6), rgb(0.8901961, 0.1019608, 0.1098039), rgb(0.9921569, 0.7490196, 0.4352941), rgb(1, 0.4980392, 0), rgb(0.7921569, 0.6980392, 0.8392157), rgb(0.4156863, 0.2392157, 0.6039216), rgb(1,0,1,0.7), rgb(0.6941176, 0.3490196, 0.1568627), rgb(0.5,0.5,0,0.8), rgb(0,0.5,0.5,0.8), rgb(0.5, 0, 0.5, 0.8), 'yellow', 'lightblue', 'pink', 'orange', 'lightgreen', 'grey30', 'grey70')
names(duplication.colours) <- paste('Dup', c(as.character(1:n.colours), '1a', '1b'), sep = '')
gene.colours <- c('black', 'purple', 'orange', 'magenta', 'brown', 'green', 'violet', 'red', 'blue')

# Create a function that can be used to plot discordant read pairs
add.diagnostics <- function(coordinates, this.col = 'red', yrange = c(0,1)){
	coordinates <- as.matrix(coordinates)
	n <- nrow(coordinates)
	co <- ncol(coordinates)
	jitter.values <- seq(yrange[1], yrange[2], length.out = n)
	points(coordinates, matrix(jitter.values, nrow = n, ncol = co), col = this.col)
	if (co == 2)
		segments(coordinates[,1], jitter.values, coordinates[,2], jitter.values, col = this.col)
}

plot.ace1 <- function(this.sample, list.of.dups = NULL, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$ace1[1], end.pos = plotting.ranges$ace1[2]){
	start.index <- which(compact.hmm$ace1[[this.sample]]$Position >= start.pos)[1]
	end.index <- tail(which(compact.hmm$ace1[[this.sample]]$Position <= end.pos), 1)
	plot(compact.hmm$ace1[[this.sample]]$Position[start.index : end.index], compact.hmm$ace1[[this.sample]]$Normalised_coverage[start.index : end.index], main = this.sample, ylim = c(0,12))
	# We highlight along the x axis where the duplications and deletions are as predicted by the FA and FM reads
	if (!is.null(list.of.dups)){
		if (list.of.dups['Dup1'])
			rect(known.cnvs$ace1$Dup1$BP$pos[1], -0.3, known.cnvs$ace1$Dup1$BP$pos[2], 0, col = rgb(1, 0, 0, 0.8), border = rgb(1, 0, 0, 0.8))
		if (list.of.dups['Dup2'])
			rect(3448100, -0.3, 3518750, 0, col = rgb(0, 1, 0, 0.8), border = rgb(0, 1, 0, 0.8))
		if (list.of.dups['Del1'])
			rect(3502000, -0.3, 3598900, 0, col = rgb(0, 0, 1, 0.8), border = rgb(0, 0, 1, 0.8))
		if (list.of.dups['Del2'])
			rect(3539450, -0.3, 3573600, 0, col = rgb(1, 0, 1, 0.8), border = rgb(1, 0, 1, 0.8))
		if (list.of.dups['Del3'])
			rect(3536000, -0.3, 3618850, 0, col = rgb(1, 1, 0, 0.8), border = rgb(1, 1, 0, 0.8))
		if (list.of.dups['Del4'])
			rect(3512500, -0.3, 3616000, 0, col = rgb(0, 1, 1, 0.8), border = rgb(0, 1, 1, 0.8))
	}
	discordant.colours <- c(FA = 'blue', SS = 'cyan', FM = 'green', XC = 'red')
	y.ranges <- list(FA = c(0,3), SS = c(3,6), FM = c(6,9), XC = c(9,12))
	for (d in diagnostics){
		if (d == 'BP'){
			add.diagnostics(diagnostic.reads$ace1[[this.sample]][[d]]$CEP$Position, this.col = 'brown', yrange = c(0,12))
			add.diagnostics(diagnostic.reads$ace1[[this.sample]][[d]]$CSP$Position, this.col = 'pink', yrange = c(0,12))
		}
		else if (d == 'XC')
			add.diagnostics(diagnostic.reads$ace1[[this.sample]][[d]]$Position, this.col = discordant.colours[d], yrange = y.ranges[[d]])
		else
			add.diagnostics(diagnostic.reads$ace1[[this.sample]][[d]], this.col = discordant.colours[d], yrange = y.ranges[[d]])
	}
	these.cnv.states <- compact.hmm$ace1[[this.sample]]$CNV[start.index : end.index]
	lines(compact.hmm$ace1[[this.sample]]$Position[start.index : end.index], these.cnv.states , col = 2, lwd = 2)
	abline(v = gene.coords$ace1[, c('start', 'end')])
	text(apply(gene.coords$ace1[, c('start', 'end')], 1, mean), 11, rownames(gene.coords$ace1), srt=90, adj = 0)
}

plot.all.ace1 <- function(list.of.samples, matrix.of.read.dups = read.based.cnvs$ace1, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$ace1[1], end.pos = plotting.ranges$ace1[2]){
	x11()
	x.midpoint <- mean(c(start.pos, end.pos))
	i <- 1
	while(1){
		if (i < 1)
			i <- length(list.of.samples)
		this.sample <- list.of.samples[i]
		if (is.null(matrix.of.read.dups))
			plot.ace1(this.sample, NULL, diagnostics, start.pos, end.pos)
		else
			plot.ace1(this.sample, matrix.of.read.dups[this.sample,], diagnostics, start.pos, end.pos)
		x <- locator(1)$x
		if (x <= x.midpoint)
			i <- i-1
		else
			i <- i+1
	}
}

plot.cyp6 <- function(this.sample, list.of.dups = NULL, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$cyp6[1], end.pos = plotting.ranges$cyp6[2]){
	start.index <- which(compact.hmm$cyp6[[this.sample]]$Position >= start.pos)[1]
	end.index <- tail(which(compact.hmm$cyp6[[this.sample]]$Position <= end.pos), 1)
	plot(compact.hmm$cyp6[[this.sample]]$Position[start.index : end.index], compact.hmm$cyp6[[this.sample]]$Normalised_coverage[start.index : end.index], main = this.sample, ylim = c(0,12))
	# We highlight along the x axis where the duplications are as predicted by the FA and SS reads
	if (!is.null(list.of.dups)){
		Dup.order <- paste('Dup', c(as.character(c(30:20, 19, 14, 15, 13:11, 18, 10:7, 17, 6:2)), '1a', '1b'), sep = '')
		list.of.dups <- list.of.dups[Dup.order]
		for (d in names(list.of.dups)[list.of.dups]){
			if (d == 'Dup15'){
				rect(known.cnvs$cyp6[[d]]$BP$pos[1], -0.6, 28555300, 0, col = duplication.colours[d], border = col)
				text(28555300, -0.4, d, cex = 1.2, adj = c(1,0))
			}
			else if (d == 'Dup23'){
				rect(28450000, -0.6, known.cnvs$cyp6[[d]]$BP$pos[2], 0, col = duplication.colours[d], border = col)
				text(known.cnvs$cyp6[[d]]$BP$pos[2], -0.4, d, cex = 1.2, adj = c(1,0))
			}
			else if (d == 'Dup27'){
				rect(28496700, -0.6, 28499200, 0, col = duplication.colours[d], border = col)
				text(28499200, -0.4, d, cex = 1.2, adj = c(1,0))
			}
			else{
				rect(known.cnvs$cyp6[[d]]$BP$pos[1], -0.6, known.cnvs$cyp6[[d]]$BP$pos[2], 0, col = duplication.colours[d], border = col)
				text(known.cnvs$cyp6[[d]]$BP$pos[2], -0.4, d, cex = 1.2, adj = c(1,0))
			}
		}
	}
	discordant.colours <- c(FA = 'blue', SS = 'cyan', FM = 'green', XC = 'red')
	y.ranges <- list(FA = c(0,3), SS = c(3,6), FM = c(6,9), XC = c(9,12))
	for (d in diagnostics){
		if (d == 'BP'){
			add.diagnostics(diagnostic.reads$cyp6[[this.sample]][[d]]$CEP$Position, this.col = 'brown', yrange = c(0,12))
			add.diagnostics(diagnostic.reads$cyp6[[this.sample]][[d]]$CSP$Position, this.col = 'pink', yrange = c(0,12))
		}
		else if (d == 'XC')
			add.diagnostics(diagnostic.reads$cyp6[[this.sample]][[d]]$Position, this.col = discordant.colours[d], yrange = y.ranges[[d]])
		else
			add.diagnostics(diagnostic.reads$cyp6[[this.sample]][[d]], this.col = discordant.colours[d], yrange = y.ranges[[d]])
	}
	these.cnv.states <- compact.hmm$cyp6[[this.sample]]$CNV[start.index : end.index]
	lines(compact.hmm$cyp6[[this.sample]]$Position[start.index : end.index], these.cnv.states , col = 2, lwd = 2)
	# Add the genes to the plot
	these.gene.colours <- gene.colours[1:nrow(gene.coords$cyp6)]
	abline(v = unlist(gene.coords$cyp6[, c('start', 'end')]), col = rep(these.gene.colours, 2))
	text(apply(gene.coords$cyp6[, c('start', 'end')], 1, mean), 9, rownames(gene.coords$cyp6), srt=90, adj = 0, col = these.gene.colours)
}

plot.all.cyp6 <- function(list.of.samples, matrix.of.read.dups = read.based.cnvs$cyp6, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$cyp6[1], end.pos = plotting.ranges$cyp6[2]){
	x11()
	x.midpoint <- mean(c(start.pos, end.pos))
	i <- 1
	while(1){
		if (i < 1)
			i <- length(list.of.samples)
		this.sample <- list.of.samples[i]
		if (is.null(matrix.of.read.dups))
			plot.cyp6(this.sample, NULL, diagnostics, start.pos, end.pos)
		else
			plot.cyp6(this.sample, matrix.of.read.dups[this.sample,], diagnostics, start.pos, end.pos)
		x <- locator(1)$x
		if (x <= x.midpoint)
			i <- i-1
		else
			i <- i+1
	}
}

plot.cyp6mz <- function(this.sample, list.of.dups = NULL, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$cyp6mz[1], end.pos = plotting.ranges$cyp6mz[2]){
	start.index <- which(compact.hmm$cyp6mz[[this.sample]]$Position >= start.pos)[1]
	end.index <- tail(which(compact.hmm$cyp6mz[[this.sample]]$Position <= end.pos), 1)
	plot(compact.hmm$cyp6mz[[this.sample]]$Position[start.index : end.index], compact.hmm$cyp6mz[[this.sample]]$Normalised_coverage[start.index : end.index], main = this.sample, ylim = c(0,12))
	# We highlight along the x axis where the duplications are as predicted by the FA and SS reads
	list.of.dups <- list.of.dups[-1]
	if (!is.null(list.of.dups)){
		these.duplication.colours <- setNames(duplication.colours[1:length(list.of.dups)], names(list.of.dups))
		for (d in names(list.of.dups)[list.of.dups]){
			if (d == 'Dupz2'){
				rect(6975100, -0.6, known.cnvs$cyp6mz$Dupz2$BP$pos[2], 0, col = these.duplication.colours[d], border = col)
				text(known.cnvs$cyp6mz$Dupz2$BP$pos[2], -0.4, d, cex = 1.2, adj = c(1,0))
			}
			else if (d == 'Dupz3'){
				rect(known.cnvs$cyp6mz$Dupz3$BP$pos[1], -0.6, 6978100, 0, col = these.duplication.colours[d], border = col)
				text(6978100, -0.4, d, cex = 1.2, adj = c(1,0))
			}
			else if (d == 'Dupz5'){
				rect(6969800, -0.6, 6976000, 0, col = these.duplication.colours[d], border = col)
				text(6976000, -0.4, d, cex = 1.2, adj = c(1,0))
			}
			else if (d == 'Dupm1'){
				rect(known.cnvs$cyp6mz$Dupm1$BP$pos[1], -0.6, 6930500, 0, col = these.duplication.colours[d], border = col)
				text(6930500, -0.4, d, cex = 1.2, adj = c(1,0))
			}
			else if (d == 'Dupm2'){
				rect(6933100, -0.6, 6935200, 0, col = these.duplication.colours[d], border = col)
				text(6935200, -0.4, d, cex = 1.2, adj = c(1,0))
			}
			else if (d == 'Dupm3'){
				rect(6929600, -0.6, 6932800, 0, col = these.duplication.colours[d], border = col)
				text(6932800, -0.4, d, cex = 1.2, adj = c(1,0))
			}
			else if (d == 'Dupmz1'){
				rect(6904000, -0.6, 6982700, 0, col = these.duplication.colours[d], border = col)
				text(6982700, -0.4, d, cex = 1.2, adj = c(1,0))
			}
			else{
				rect(known.cnvs$cyp6mz[[d]]$BP$pos[1], -0.6, known.cnvs$cyp6mz[[d]]$BP$pos[2], 0, col = these.duplication.colours[d], border = col)
				text(known.cnvs$cyp6mz[[d]]$BP$pos[2], -0.4, d, cex = 1.2, adj = c(1,0))
			}
		}
	}
	discordant.colours <- c(FA = 'blue', SS = 'cyan', FM = 'green', XC = 'red')
	y.ranges <- list(FA = c(0,3), SS = c(3,6), FM = c(6,9), XC = c(9,12))
	for (d in diagnostics){
		if (d == 'BP'){
			add.diagnostics(diagnostic.reads$cyp6mz[[this.sample]][[d]]$CEP$Position, this.col = 'brown', yrange = c(0,12))
			add.diagnostics(diagnostic.reads$cyp6mz[[this.sample]][[d]]$CSP$Position, this.col = 'pink', yrange = c(0,12))
		}
		else if (d == 'XC')
			add.diagnostics(diagnostic.reads$cyp6mz[[this.sample]][[d]]$Position, this.col = discordant.colours[d], yrange = y.ranges[[d]])
		else
			add.diagnostics(diagnostic.reads$cyp6mz[[this.sample]][[d]], this.col = discordant.colours[d], yrange = y.ranges[[d]])
	}
	these.cnv.states <- compact.hmm$cyp6mz[[this.sample]]$CNV[start.index : end.index]
	lines(compact.hmm$cyp6mz[[this.sample]]$Position[start.index : end.index], these.cnv.states , col = 2, lwd = 2)
	# Add the genes to the plot
	these.gene.colours <- gene.colours[1:nrow(gene.coords$cyp6mz)]
	abline(v = unlist(gene.coords$cyp6mz[, c('start', 'end')]), col = rep(these.gene.colours, 2))
	text(apply(gene.coords$cyp6mz[, c('start', 'end')], 1, mean), 9, rownames(gene.coords$cyp6mz), srt=90, adj = 0, col = these.gene.colours)
}

plot.all.cyp6mz <- function(list.of.samples, matrix.of.read.dups = read.based.cnvs$cyp6mz, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$cyp6mz[1], end.pos = plotting.ranges$cyp6mz[2]){
	x11()
	x.midpoint <- mean(c(start.pos, end.pos))
	i <- 1
	while(1){
		if (i < 1)
			i <- length(list.of.samples)
		this.sample <- list.of.samples[i]
		if (is.null(matrix.of.read.dups))
			plot.cyp6mz(this.sample, NULL, diagnostics, start.pos, end.pos)
		else
			plot.cyp6mz(this.sample, matrix.of.read.dups[this.sample,], diagnostics, start.pos, end.pos)
		x <- locator(1)$x
		if (x <= x.midpoint)
			i <- i-1
		else
			i <- i+1
	}
}

plot.gste <- function(this.sample, list.of.dups = NULL, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$gste[1], end.pos = plotting.ranges$gste[2]){
	start.index <- which(compact.hmm$gste[[this.sample]]$Position >= start.pos)[1]
	end.index <- tail(which(compact.hmm$gste[[this.sample]]$Position <= end.pos), 1)
	plot(compact.hmm$gste[[this.sample]]$Position[start.index : end.index], compact.hmm$gste[[this.sample]]$Normalised_coverage[start.index : end.index], main = this.sample, ylim = c(0,12))
	# We highlight along the x axis where the duplications are as predicted by the FA and SS reads
	if (!is.null(list.of.dups)){
		# The order in which we plot Dups is important for clarity (we want smaller ones on top of larger ones)
		list.of.dups <- rev(list.of.dups[-1])
		for (d in names(list.of.dups)[list.of.dups]){
			if (d == 'Dup13'){
				rect(known.cnvs$gste[[d]]$BP$pos[1], -0.6, 28599287, 0, col = duplication.colours[d], border = col)
				text(28599287, -0.4, d, cex = 1.2, adj = c(1,0))
			}
			else{
				rect(known.cnvs$gste[[d]]$BP$pos[1], -0.6, known.cnvs$gste[[d]]$BP$pos[2], 0, col = duplication.colours[d], border = col)
				text(known.cnvs$gste[[d]]$BP$pos[2], -0.4, d, cex = 1.2, adj = c(1,0))
			}
		}
	}
	discordant.colours <- c(FA = 'blue', SS = 'cyan', FM = 'green', XC = 'red')
	y.ranges <- list(FA = c(0,3), SS = c(3,6), FM = c(6,9), XC = c(9,12))
	for (d in diagnostics){
		if (d == 'BP'){
			add.diagnostics(diagnostic.reads$gste[[this.sample]][[d]]$CEP$Position, this.col = 'brown', yrange = c(0,12))
			add.diagnostics(diagnostic.reads$gste[[this.sample]][[d]]$CSP$Position, this.col = 'pink', yrange = c(0,12))
		}
		else if (d == 'XC')
			add.diagnostics(diagnostic.reads$gste[[this.sample]][[d]]$Position, this.col = discordant.colours[d], yrange = y.ranges[[d]])
		else
			add.diagnostics(diagnostic.reads$gste[[this.sample]][[d]], this.col = discordant.colours[d], yrange = y.ranges[[d]])
	}
	these.cnv.states <- compact.hmm$gste[[this.sample]]$CNV[start.index : end.index]
	lines(compact.hmm$gste[[this.sample]]$Position[start.index : end.index], these.cnv.states , col = 2, lwd = 2)
	# Add the genes to the plot
	these.gene.colours <- gene.colours[1:nrow(gene.coords$gste)]
	abline(v = unlist(gene.coords$gste[, c('start', 'end')]), col = rep(these.gene.colours, 2))
	text(apply(gene.coords$gste[, c('start', 'end')], 1, mean), 9, rownames(gene.coords$gste), srt=90, adj = 0, col = these.gene.colours)
}

plot.all.gste <- function(list.of.samples, matrix.of.read.dups = read.based.cnvs$gste, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$gste[1], end.pos = plotting.ranges$gste[2]){
	x11()
	x.midpoint <- mean(c(start.pos, end.pos))
	i <- 1
	while(1){
		if (i < 1)
			i <- length(list.of.samples)
		this.sample <- list.of.samples[i]
		if (is.null(matrix.of.read.dups))
			plot.gste(this.sample, NULL, diagnostics, start.pos, end.pos)
		else
			plot.gste(this.sample, matrix.of.read.dups[this.sample,], diagnostics, start.pos, end.pos)
		x <- locator(1)$x
		if (x <= x.midpoint)
			i <- i-1
		else
			i <- i+1
	}
}

plot.cyp9k1 <- function(this.sample, list.of.dups = NULL, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$cyp9k1[1], end.pos = plotting.ranges$cyp9k1[2]){
	start.index <- which(compact.hmm$cyp9k1[[this.sample]]$Position >= start.pos)[1]
	end.index <- tail(which(compact.hmm$cyp9k1[[this.sample]]$Position <= end.pos), 1)
	plot(compact.hmm$cyp9k1[[this.sample]]$Position[start.index : end.index], compact.hmm$cyp9k1[[this.sample]]$Normalised_coverage[start.index : end.index], main = this.sample, ylim = c(0,12))
	# We highlight along the x axis where the duplications are as predicted by the FA and SS reads
	if (!is.null(list.of.dups)){
		# The order in which we plot Dups is important for clarity (we want smaller ones on top of larger ones)
		Dup.order <- paste('Dup', as.character(c(28:1, 17)), sep = '')
		list.of.dups <- list.of.dups[Dup.order]
		for (d in names(list.of.dups)[list.of.dups]){
			if (d == 'Dup3'){
				rect(15240464, -0.6, known.cnvs$cyp9k1[[d]]$BP$pos[2], 0, col = duplication.colours[d], border = col)
				text(known.cnvs$cyp9k1[[d]]$BP$pos[2], -0.4, d, cex = 1.2, adj = c(1,0))
			}
			else if (d == 'Dup7'){
				rect(15238000, -0.6, 15246300,  0, col = duplication.colours[d], border = col)
				text(15246300, -0.4, d, cex = 1.2, adj = c(1,0))
			}
			else if (d == 'Dup14'){
				rect(15233807, -0.6, known.cnvs$cyp9k1[[d]]$BP$pos[2], 0, col = duplication.colours[d], border = col)
				text(known.cnvs$cyp9k1[[d]]$BP$pos[2], -0.4, d, cex = 1.2, adj = c(1,0))
			}
			else if (d == 'Dup15'){
				rect(15233807, -0.6, known.cnvs$cyp9k1[[d]]$BP$pos[2], 0, col = duplication.colours[d], border = col)
				text(known.cnvs$cyp9k1[[d]]$BP$pos[2], -0.4, d, cex = 1.2, adj = c(1,0))
			}
			else if (d == 'Dup16'){
				rect(15222810, -0.6, known.cnvs$cyp9k1[[d]]$BP$pos[2], 0, col = duplication.colours[d], border = col)
				text(known.cnvs$cyp9k1[[d]]$BP$pos[2], -0.4, d, cex = 1.2, adj = c(1,0))
			}
			else if (d == 'Dup18'){
				rect(15236100, -0.6, known.cnvs$cyp9k1[[d]]$BP$pos[2], 0, col = duplication.colours[d], border = col)
				text(known.cnvs$cyp9k1[[d]]$BP$pos[2], -0.4, d, cex = 1.2, adj = c(1,0))
			}
			else if (d == 'Dup24'){
				rect(15238550, -0.6, 15255100, 0, col = duplication.colours[d], border = col)
				text(15255100, -0.4, d, cex = 1.2, adj = c(1,0))
			}
			else if (d == 'Dup25'){
				rect(15223500, -0.6, 15246650, 0, col = duplication.colours[d], border = col)
				text(15246650, -0.4, d, cex = 1.2, adj = c(1,0))
			}
			else if (d == 'Dup26'){
				rect(15222700, -0.6, 15248050, 0, col = duplication.colours[d], border = col)
				text(15248050, -0.4, d, cex = 1.2, adj = c(1,0))
			}
			else if (d == 'Dup27'){
				rect(known.cnvs$cyp9k1[[d]]$BP$pos[1], -0.6, 15248650, 0, col = duplication.colours[d], border = col)
				text(15248650, -0.4, d, cex = 1.2, adj = c(1,0))
			}
			else if (d == 'Dup28'){
				rect(15235000, -0.6, known.cnvs$cyp9k1[[d]]$BP$pos[2], 0, col = duplication.colours[d], border = col)
				text(known.cnvs$cyp9k1[[d]]$BP$pos[2], -0.4, d, cex = 1.2, adj = c(1,0))
			}
			else{
				rect(known.cnvs$cyp9k1[[d]]$BP$pos[1], -0.6, known.cnvs$cyp9k1[[d]]$BP$pos[2], 0, col = duplication.colours[d], border = col)
				text(known.cnvs$cyp9k1[[d]]$BP$pos[2], -0.4, d, cex = 1.2, adj = c(1,0))
			}
		}
	}
	discordant.colours <- c(FA = 'blue', SS = 'cyan', FM = 'green', XC = 'red')
	y.ranges <- list(FA = c(0,3), SS = c(3,6), FM = c(6,9), XC = c(9,12))
	for (d in diagnostics){
		if (d == 'BP'){
			add.diagnostics(diagnostic.reads$cyp9k1[[this.sample]][[d]]$CEP$Position, this.col = 'brown', yrange = c(0,12))
			add.diagnostics(diagnostic.reads$cyp9k1[[this.sample]][[d]]$CSP$Position, this.col = 'pink', yrange = c(0,12))
		}
		else if (d == 'XC')
			add.diagnostics(diagnostic.reads$cyp9k1[[this.sample]][[d]]$Position, this.col = discordant.colours[d], yrange = y.ranges[[d]])
		else
			add.diagnostics(diagnostic.reads$cyp9k1[[this.sample]][[d]], this.col = discordant.colours[d], yrange = y.ranges[[d]])
	}
	these.cnv.states <- compact.hmm$cyp9k1[[this.sample]]$CNV[start.index : end.index]
	lines(compact.hmm$cyp9k1[[this.sample]]$Position[start.index : end.index], these.cnv.states , col = 2, lwd = 2)
	# Add the genes to the plot
	abline(v = c(gene.coords$cyp9k1[, c('start', 'end')]))
	text(apply(gene.coords$cyp9k1[, c('start', 'end')], 1, mean), 11, rownames(gene.coords$cyp9k1), srt=90, adj = 0)
}

plot.all.cyp9k1 <- function(list.of.samples, matrix.of.read.dups = read.based.cnvs$cyp9k1, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$cyp9k1[1], end.pos = plotting.ranges$cyp9k1[2]){
	x11()
	x.midpoint <- mean(c(start.pos, end.pos))
	i <- 1
	while(1){
		if (i < 1)
			i <- length(list.of.samples)
		this.sample <- list.of.samples[i]
		if (is.null(matrix.of.read.dups))
			plot.cyp9k1(this.sample, NULL, diagnostics, start.pos, end.pos)
		else
			plot.cyp9k1(this.sample, matrix.of.read.dups[this.sample,], diagnostics, start.pos, end.pos)
		x <- locator(1)$x
		if (x <= x.midpoint)
			i <- i-1
		else
			i <- i+1
	}
}

simple.plot <- function(sample.names, gene.cluster, smoothing = 5, maxy = 12, ...){
	these.coords <- plotting.ranges[[gene.cluster]]
	plot(c(min(these.coords), max(these.coords)), c(0,maxy), type = 'n', bty = 'n', xaxt = 'n', xlab = '', ylab = 'Normalised coverage')
	axis(1, lwd = 0, mgp = c(0,0.25,2))
	mtext(paste('Position on chromosome', unique(gene.coords[[gene.cluster]]$Chrom)), 1, 1.75)
	these.gene.colours <- gene.colours[1:nrow(gene.coords[[gene.cluster]])]
	abline(v = unlist(gene.coords[[gene.cluster]][, c('start', 'end')]), col = rep(these.gene.colours, 2))
	text(apply(gene.coords[[gene.cluster]][, c('start', 'end')], 1, mean), 9, rownames(gene.coords[[gene.cluster]]), srt=90, adj = 0, col = these.gene.colours)
	smoothed.line <- function(s){
		these.counts <- subset(compact.hmm[[gene.cluster]][[s]], Position >= these.coords[1] & Position <= these.coords[2])
		# If needed, create a table of smoothed data here
		if (smoothing > 1){
			if (smoothing > nrow(these.counts))
				stop('Fail. Smoothing window size is larger than the number of points to smooth over.')
			these.counts <- t(sapply(1:(nrow(these.counts) - smoothing + 1), function(j) apply(these.counts[j:(j + smoothing - 1), c('Position', 'Normalised_coverage')], 2, mean)))
		}
		lines(these.counts[, 'Position'], these.counts[, 'Normalised_coverage'])
	}
	empty <- sapply(sample.names, smoothed.line)
}

