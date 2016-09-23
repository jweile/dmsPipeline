#######################################
# Codon preference vs Complementation #
#######################################

options(stringsAsFactors=FALSE)

source("lib/resultfile.R")
source("lib/libyogitools.R")
source("lib/liblogging.R")
source("lib/cliargs.R")

library("beeswarm")


#Get output directory
outdir <- getArg("outdir",default="workspace/test/")

#Set resultfile
html <- new.resultfile(paste0(outdir,"results.html"))
html$section("Codon preference analysis")

#Init logger
logger <- new.logger(paste0(outdir,"codonPreference.log"))

logger$info("Reading input data")

codonPref <- read.delim("res/codonPref.tsv")
rownames(codonPref) <- codonPref$Codon

codonPref$pref <- sapply(1:nrow(codonPref), function(i) {
	aa <- codonPref[i,"AmAcid"]
	n <- sum(codonPref[,"AmAcid"]==aa)
	codonPref[i,"Fraction"]-(1/n)
})

orf.start <- 76
orf.end <- 549
amplicon <- scan("res/ube2i-comp.fa",what="character")[[2]]
orf <- substr(amplicon,orf.start,orf.end)
codon.starts <- seq(1,nchar(orf),3)
codons <- sapply(codon.starts,function(i)substr(orf,i,i+2))


ube2i.pref <- data.frame(pos=1:length(codons),bcodon=codons, pref=codonPref[codons,"pref"])

screen <- read.csv(paste0(outdir,"compl_joint_results_UBE2I.csv"))
screen$pos <- as.integer(substr(screen$mut,2,nchar(screen$mut)-1))

posAverages <- with(screen, tapply(score,pos,mean,na.rm=TRUE))

ube2i.pref$posAverage <- posAverages[ube2i.pref$pos]

ube2i.pref$category <- factor(sapply(ube2i.pref$posAverage,function(x){
	if (x < .5) "s < 0.5"
	else if (x < 1.2) "0.5 < s < 1.2"
	else "s > 1.2"
}),levels=c("s < 0.5","0.5 < s < 1.2","s > 1.2"))

logger$info("Drawing plot")

html$figure(function(){
	with(ube2i.pref,beeswarm(jitter(pref) ~ category,
		col="steelblue3",pch=20,method="swarm",
		ylab="Yeast codon preference",
		xlab="Complementation growth score"
	))
	with(ube2i.pref,bxplot(pref ~ category,add=TRUE))
},paste0(outdir,"codonPref"),5,5)

html$shutdown()

logger$info("Done!")
