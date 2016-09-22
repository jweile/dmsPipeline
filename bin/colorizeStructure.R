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

logger$info("Calculating position-wise statistics")

data.stats <- do.call(rbind,tapply(data$joint.score,pos,summary))
colnames(data.stats) <- c("min","q25","med","mean","q75","max")

logger$info("Colorizing structure")

invisible(lapply(c("min","med","max"),function(stat) {
	outfile <- paste0(outdir,"compl_",geneName,"_pymol_colors_",stat,".txt")
	pycol <- new.pymol.colorizer(outfile)
	pycol$define.colors()
	pycol$colorize(cbind(
		index=as.integer(rownames(data.stats)),
		fitness=data.stats[,stat]
	))
	pycol$close()
	html$link.data(outfile)
}))

logger$info("Done")
