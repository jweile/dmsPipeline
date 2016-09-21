#####################################
# Compare Y2H to compl data to find #
# interface positions               #
#####################################

source("lib/resultfile.R")
source("lib/liblogging.R")
source("lib/libyogitools.R")
source("lib/cliargs.R")
source("lib/topoScatter.R")

library("hash")
options(stringsAsFactors=FALSE)


#get output directory
outdir <- getArg("outdir",default="workspace/test/")

#Initialize logger
logger <- new.logger(paste0(outdir,"findInterfaces.log"))

#Set resultfile
html <- new.resultfile(paste0(outdir,"results.html"))
html$section("UBE2I Y2H BarSEQ analysis")


logger$info("Reading input")

y2h <- read.csv(paste0(outdir,"y2h_scores_perMut.csv"))
compl <- read.csv(paste0(outdir,"imputed_regularized_UBE2I_scores.csv"))
rownames(compl) <- compl$mut

#remove multi-mutants
data <- y2h[regexpr(",",y2h$mut) < 0,]

data$compl.score <- compl[data$mut,"joint.score"]
data$compl.sd <- compl[data$mut,"joint.sd"]

interactors <- c("SATB1","SUMO1","ZBED1","RNF40")
candidates <- NULL

logger$info("Comparing Y2H to Complementation scores")

html$subsection("Complementation vs Interaction Scores")
html$figure(function() {
	op <- par(mfrow=c(2,2))
	candidates <<- lapply(interactors,function(ia) {
		iasd <- paste0(ia,".sd")
		iav <- paste0(ia,".score")
		is <- which(!is.na(data[,iasd]) & 
			data[,iasd] > 0 & 
			data[,iasd] < 0.5 & 
			data$compl.sd < 0.3
		)
		plot(NA,type="n",
			main=paste("Complementation vs",ia),
			xlab="Complementation score",
			ylab=paste(ia,"interaction score"),
			xlim=c(-0.5,2),
			ylim=c(-1,2)
		)
		x <- data[is,c("compl.score")]
		y <- data[is,c(iav)]
		xsd <- data[is,c("compl.sd")]
		ysd <- data[is,c(iasd)]
		arrows(x-xsd/2,y,x+xsd/2,y,length=0.01,angle=90,code=3)
		arrows(x,y-ysd/2,x,y+ysd/2,length=0.01,angle=90,code=3)
		abline(h=0:1,v=0:1,col=c("firebrick3","chartreuse3"))

		m <- data[is,"mut"]
		js <- which(x > 0.5 &  y < 0.5 & y < x-0.5)
		cand.j <- data.frame(mut=m[js],compl=x[js],ia=y[js],type="interfacial")
		ks <- which(x < 0.5 &  y > 0.5 & y > x+0.5)
		cand.k <- data.frame(mut=m[ks],compl=x[ks],ia=y[ks],type="deadfold")
		rbind(cand.j,cand.k)
	})
	par(op)
},paste0(outdir,"complVinteraction"),10,10)

names(candidates) <- interactors

outfile <- paste0(outdir,"interface_candidates.txt")
con <- file(outfile,open="w")
invisible(lapply(interactors,function(ia) {
	writeLines(c(ia,"======="),con)
	write.table(format(candidates[[ia]],digits=2),con,quote=FALSE,sep="\t",row.names=FALSE)
	writeLines("\n\n",con)
}))
close(con)

html$subsection("Interface candidates")
html$link.data(outfile)

html$shutdown()

logger$info("Done!")

