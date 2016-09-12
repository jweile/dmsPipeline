
#####################################################
# Plot fitness scores against surface accessibility #
# and for core/surface/interface residues           #
#####################################################

source("lib/resultfile.R")
source("lib/liblogging.R")
source("lib/cliargs.R")

# library("beeswarm")
options(stringsAsFactors=FALSE)

outdir <- getArg("outdir",default="workspace/test/")

#Set resultfile
html <- new.resultfile(paste0(outdir,"results.html"))
html$section("Surface accessibility and conservation analysis")

#Init logger
logger <- new.logger(paste0(outdir,"accessibility.log"))

##############
# LOAD INPUT #
##############
logger$info("Reading input")


struc <- read.csv("res/UBE2I_structural_features.csv")
rownames(struc) <- struc$pos

data <- read.csv(paste0(outdir,"imputed_regularized_UBE2I_scores.csv"))
singles <- data[regexpr(",",data$mut) < 1,]
singles <- singles[regexpr("\\w\\d+\\w",singles$mut) > 0,]

##############
# CALCULATE #
##############

singles.pos <- as.numeric(substr(singles$mut,2,nchar(singles$mut)-1))

logger$info("Looking up solvent accessibility for each mutant residue")

singles.acc <- struc[as.character(singles.pos),"r.acc"]

logger$info("Looking up maximal burial for each mutant residue")

singles.maxBurial <- apply(
	struc[as.character(singles.pos),c("b.sumo","b.sumo.nc","b.e1","b.rangap","b.ranbp2","b.homodimer")],
	1,max,na.rm=TRUE
)
singles.maxBurial[is.infinite(singles.maxBurial)] <- NA

logger$info("Deviding mutations into accessibility classes")

acc.groups <- list(
	core=singles[singles.acc < 0.2,"joint.score"],
	interface=singles[singles.maxBurial > 0.4,"joint.score"],
	surface=singles[singles.acc > 0.6 & singles.maxBurial < 0.2,"joint.score"]
)

logger$info("Deviding mutations into conservation classes")

ube2i.cons <- as.integer(scan("res/UBE2I_cons.txt",what="integer",sep=","))
singles.cons <- ube2i.cons[singles.pos]

cons.groups <- list(
	low=singles[singles.cons < 4,"joint.score"],
	medium=singles[singles.cons >= 4 & singles.cons < 7,"joint.score"],
	high=singles[singles.cons >= 7,"joint.score"]
)

########################################
# DRAW PLOT AND PERFORM WILCOXON TESTS #
########################################

logger$info("Drawing plot")

# pdf(paste0(outdir,"accessibilityAndConservation.pdf"),7.5,4)
html$figure(function(){
	layout(cbind(1,2))
	xs <- boxplot(
		acc.groups,
		col=c("firebrick3","steelblue3","chartreuse3"),
		ylab="fitness score",xlab="Residue accessibility",
		ylim=c(-0.5,4)
	)
	coreSurf <- with(acc.groups,wilcox.test(core,surface,alternative="less"))
	if (coreSurf$p.value < 0.05) {
		lines(c(1,1,3,3),c(3.4,3.5,3.5,3.4))
		text(2,3.6,"*",cex=1.4)
	}
	ifSurf <- with(acc.groups,wilcox.test(interface,surface,alternative="less"))
	if (ifSurf$p.value < 0.05) {
		lines(c(2,2,3,3),c(3,3.1,3.1,3))
		text(2.5,3.2,"*",cex=1.4)
	}

	boxplot(
		cons.groups,
		col=c("firebrick3","orange","gold1"),
		ylab="fitness score",xlab="Evolutionary conservation",
		ylim=c(-0.5,3.5)
	)
	lo.mid <- with(cons.groups,wilcox.test(low,medium,alternative="greater"))
	if (lo.mid$p.value < 0.05) {
		y <- 3.15
		lines(c(1,1,2,2),c(y-.1,y,y,y-.1))
		text(1.5,y+.1,"*",cex=1.4)
	}

	mid.hi <- with(cons.groups,wilcox.test(medium,high,alternative="greater"))
	if (mid.hi$p.value < 0.05) {
		y <- 3
		lines(c(2,2,3,3),c(y-.1,y,y,y-.1))
		text(2.5,y+.1,"*",cex=1.4)
	}
},paste0(outdir,"accessibilityAndConservation"),7.5,4)
# invisible(dev.off())

html$shutdown()

logger$info("Done.")
