#################################
# Colorize structures based on  #
# complementation data          #
#################################

source("lib/resultfile.R")
source("lib/liblogging.R")
source("lib/cliargs.R")
source("lib/pymolCol.R")

library("hash")
options(stringsAsFactors=FALSE)

#get output directory
outdir <- getArg("outdir",default="workspace/test/")
infile <- getArg("infile",default="workspace/test/imputed_regularized_UBE2I_scores.csv")
geneName <- getArg("geneName",default="UBE2I")

chain <- getArg("chain",default="A")

bendParam <- as.numeric(getArg("bend",default=0))

#Initialize logger
logger <- new.logger(paste0(outdir,"colorizeStructure_",geneName,".log"))

#Set resultfile
html <- new.resultfile(paste0(outdir,"results.html"))
html$section(paste("Colorizing Structure of",geneName))

logger$info("Filtering input")

data <- read.csv(infile)
data <- data[regexpr("_",data$mut) < 1,]
data <- data[with(data,substr(mut,1,1) != substr(mut,nchar(mut),nchar(mut))),]

pos <- as.integer(substr(data$mut,2,nchar(data$mut)-1))


#carefully bend scores into a given direction
bend <- function(x,a=0) {
	if (is.na(x)) return(NA)
	if (x < 0) x <- 0
	if (x > 1) x <- 1/x
	if (a == 0) return(x)
	if (abs(a) > 0.5) stop("parameter a must be between -0.5 and 0.5")
	q <- (1 - 2*a) / (4*a)
	t <- ifelse(a >= 0,1,-1) * sqrt(x/(2*a) + q^2) - q
	(1+2*a)*t - 2*a*t^2
}

if (bendParam != 0) {
	logger$info("Applying bending transformation")
	data$joint.score <- sapply(data$joint.score,bend,a=bendParam)
}

logger$info("Calculating position-wise statistics")

data.stats <- do.call(rbind,tapply(data$joint.score,pos,summary))
colnames(data.stats) <- c("min","q25","med","mean","q75","max")

logger$info("Colorizing structure")

invisible(lapply(c("min","med","max"),function(stat) {
	outfile <- paste0(outdir,"compl_",geneName,"_pymol_colors_",stat,".txt")
	pycol <- new.pymol.colorizer(outfile,chain=chain)
	pycol$define.colors()
	pycol$colorize(cbind(
		index=as.integer(rownames(data.stats)),
		fitness=data.stats[,stat]
	))
	pycol$close()
	html$link.data(outfile)
}))

html$shutdown()

logger$info("Done.")
