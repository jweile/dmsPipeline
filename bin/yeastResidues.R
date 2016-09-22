
source("lib/resultfile.R")
source("lib/liblogging.R")
source("lib/cliargs.R")

library("beeswarm")

options(stringsAsFactors=FALSE)

#get output directory
outdir <- getArg("outdir",default="workspace/test/")

#Initialize logger
logger <- new.logger(paste0(outdir,"yeastResidues.log"))

#Set resultfile
html <- new.resultfile(paste0(outdir,"results.html"))
html$section("Yeast residue bias")

logger$info("Reading input")

compl <- read.csv(paste0(outdir,"imputed_regularized_UBE2I_scores.csv"))
compl$pos <- as.integer(substr(compl$mut,2,nchar(compl$mut)-1))
rownames(compl) <- compl$mut

foo <- scan("res/UBE2I_alignment.fa",what="character",sep="\n")
seqs <- sapply(0:12,function(i) paste(foo[i*4 + (2:4)],collapse=""))
names(seqs) <- sapply(strsplit(foo[0:12*4+1],"/"),function(x)x[[1]])

to.chars <- function(x) sapply(1:nchar(x),function(i)substr(x,i,i))

human.seq <- to.chars(seqs[[">Human"]])
yeast.seq <- to.chars(seqs[[">BakersYeast"]])
mold.seq <- to.chars(seqs[[">SlimeMold"]])
fly.seq <- to.chars(seqs[[">Fly"]])

positions <- unique(compl$pos)

yeast.scores <- na.omit(sapply(positions,function(i) {
	h <- human.seq[[i]]
	y <- yeast.seq[[i]]
	if (h != y) {
		m <- paste0(h,i,y)
		s <- compl[m,"screen.score"]
		# if (!is.na(s) && s > 1.5) cat(m,s,"\n")
		s
	} else NA
}))

mold.scores <- na.omit(sapply(positions,function(i) {
	h <- human.seq[[i]]
	y <- mold.seq[[i]]
	if (h != y) {
		m <- paste0(h,i,y)
		s <- compl[m,"screen.score"]
		# if (!is.na(s) && s > 1.5) cat(m,s,"\n")
		s
	} else NA
}))

fly.scores <- na.omit(sapply(positions,function(i) {
	h <- human.seq[[i]]
	y <- fly.seq[[i]]
	if (h != y) {
		m <- paste0(h,i,y)
		s <- compl[m,"screen.score"]
		# if (!is.na(s) && s > 1.5) cat(m,s,"\n")
		s
	} else NA
}))

logger$info("Drawing plot")

# pdf("yeastAAs.pdf",5.5,5)
html$figure(function(){
	beedata <- list(Yeast=yeast.scores,SlimeMold=mold.scores,Fly=fly.scores)
	beeswarm(beedata,
		col=c("firebrick3","chartreuse3","steelblue3"),
		labels=c("S. cerevisiae\n(Yeast)","D. discoideum\n(Slime Mold)","D. melanogaster\n(Fly)")
		,pch=16,ylab="Complementation fitness"
	)
	bxplot(beedata,add=TRUE)
	abline(h=0:1,col=c("firebrick4","darkolivegreen4"),lty="dashed")
},paste0(outdir,"yeastAAs"),5.5,5)
# dev.off()

html$shutdown()

logger$info("Done")
