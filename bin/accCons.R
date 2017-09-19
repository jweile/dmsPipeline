
#####################################################
# Plot fitness scores against surface accessibility #
# and for core/surface/interface residues           #
#####################################################

source("lib/resultfile.R")
source("lib/liblogging.R")
source("lib/cliargs.R")

library("Biostrings")
library("beeswarm")

data(BLOSUM62)

# library("beeswarm")
options(stringsAsFactors=FALSE)

outdir <- getArg("outdir",default="workspace/test/")

#Set resultfile
html <- new.resultfile(paste0(outdir,"results.html"))
html$section("Accessibility, conservation, and predictor comparison")

#Init logger
logger <- new.logger(paste0(outdir,"accessibility.log"))

##############
# LOAD INPUT #
##############
logger$info("Reading input")


struc <- read.csv("res/UBE2I_structural_features.csv")
rownames(struc) <- struc$pos

pp2 <- read.delim("res/UBE2I_polyphen2.txt",comment.char="#")
pp2$mut <- with(pp2,paste0(aa1,pos,aa2))
rownames(pp2) <- pp2$mut

provean <- read.delim("res/UBE2I_provean.tsv")
provean <- provean[!is.na(provean$POSITION),]
provean$mut <- with(provean,paste0(RESIDUE_REF,POSITION,RESIDUE_ALT))
rownames(provean) <- provean$mut

scinotExpr <- function(x, digits=2) {
    sign <- ""
    if (x < 0) {
        sign <- "-"
        x <- -x
    }
    exponent <- floor(log10(x))
    if (exponent) {
        xx <- round(x / 10^exponent, digits=digits)
        e <- paste(" %*% 10^", as.integer(exponent), sep="")
    } else {
        xx <- round(x, digits=digits)
        e <- ""
    }
    parse(text=paste("P == ",sign, xx, e, sep=""))
}

drawStar <- function(pval,x1,x2,y) {
	if (pval < 0.05) {
		lines(c(x1,x1,x2,x2),c(y-.1,y,y,y-.1))
		# text(mean(c(x1,x2)),y+.1,"*",cex=1.4)
		text(mean(c(x1,x2)),y+.1,scinotExpr(pval),cex=0.8)
	}
}

##############
# CALCULATE #
##############

calculateAndPlot <- function(singles,datalabel) {

	singles.pos <- as.numeric(substr(singles$mut,2,nchar(singles$mut)-1))

	logger$info("Looking up solvent accessibility for each mutant residue")

	singles.acc <- struc[as.character(singles.pos),"r.acc"]

	logger$info("Looking up maximal burial for each mutant residue")

	singles.maxBurial <- apply(
		struc[as.character(singles.pos),c("b.sumo","b.sumo.nc","b.e1","b.rangap","b.ranbp2","b.homodimer")],
		1,max,na.rm=TRUE
	)
	singles.maxBurial[is.infinite(singles.maxBurial)] <- NA

	logger$info("Dividing mutations into accessibility classes")

	acc.groups <- list(
		core=singles[singles.acc < 0.2,"score"],
		interface=singles[singles.maxBurial > 0.4,"score"],
		surface=singles[singles.acc > 0.6 & singles.maxBurial < 0.2,"score"]
	)

	iBurials <- struc[as.character(singles.pos),c("b.sumo","b.sumo.nc","b.e1","b.rangap","b.ranbp2","b.homodimer")]
	interface.groups <- c(lapply(colnames(iBurials),function(ifID) {
		singles[iBurials[,ifID] > 0.4,"score"]
	}),list(acc.groups$surface))
	names(interface.groups) <- c("SUMO","SUMO-nc","E1","Substrate","E3","Homodimer","None")

	logger$info("Dividing mutations into conservation classes")

	ube2i.cons <- as.integer(scan("res/UBE2I_cons.txt",what="integer",sep=","))
	singles.cons <- ube2i.cons[singles.pos]

	cons.groups <- list(
		low=singles[singles.cons < 4,"score"],
		medium=singles[singles.cons >= 4 & singles.cons < 7,"score"],
		high=singles[singles.cons >= 7,"score"]
	)

	logger$info("Looking up BLOSUM62 scores")

	singles.from <- substr(singles$mut,1,1)
	singles.to <- substr(singles$mut,nchar(singles$mut),nchar(singles$mut))
	singles.blosum <- mapply(function(from,to) BLOSUM62[from,to],from=singles.from,to=singles.to)

	logger$info("Dividing mutations into BLOSUM classes")

	blosum.groups <- list(
		low=singles[singles.blosum < 0,"score"],
		high=singles[singles.blosum >= 0,"score"]
	)

	logger$info("Dividing mutations into Polyphen2 classes")

	singles.pp2 <- pp2[singles$mut,"prediction"]
	pp2.groups <- tapply(singles$score,singles.pp2,function(x)x)

	logger$info("Dividing mutations into PROVEAN classes")

	singles.prov <- provean[singles$mut,"PREDICTION..cutoff..2.5."]
	prov.groups <- tapply(singles$score,singles.prov,function(x)x)

	logger$info("Dividing mutations into SIFT classes")

	singles.sift <- provean[singles$mut,"PREDICTION..cutoff.0.05."]
	sift.groups <- tapply(singles$score,singles.sift,function(x)x)

	########################################
	# DRAW PLOT AND PERFORM WILCOXON TESTS #
	########################################

	logger$info("Drawing plot")

	html$subsection(paste(datalabel,"data"))

	html$figure(function(){

		beeswarm(
			interface.groups,
			cex=0.6,pch=20,
			col=c("firebrick3","firebrick4","steelblue3","chartreuse3","gray60","gold2","gray40"),
			ylim=c(-0.5,2.7),
			ylab="mutant fitness score"
		)
		abline(h=0:1,col=c("firebrick2","darkolivegreen3"))
		bxplot(interface.groups,add=TRUE,col="gray30")

		cVnc <- with(interface.groups,wilcox.test(SUMO,`SUMO-nc`,alternative="greater"))
		drawStar(cVnc$p.value,1,2,2.2)
		ncVe1 <- with(interface.groups,wilcox.test(`SUMO-nc`,E1,alternative="less"))
		drawStar(ncVe1$p.value,2,3,2.2)#n.s.
		sumoVe1 <- with(interface.groups,wilcox.test(SUMO,E1,alternative="greater"))
		drawStar(sumoVe1$p.value,1,2.98,2.4)
		e1VSub <- with(interface.groups,wilcox.test(Substrate,E1,alternative="greater"))
		drawStar(e1VSub$p.value,3,4,2.2)
		noneVhomo <- with(interface.groups,wilcox.test(None,Homodimer,alternative="greater"))
		drawStar(noneVhomo$p.value,6.02,7,2.4)
		homoVe1 <- with(interface.groups,wilcox.test(Homodimer,E1,alternative="greater"))
		drawStar(homoVe1$p.value,3.02,5.98,2.4)
		homoVe3 <- with(interface.groups,wilcox.test(Homodimer,E3,alternative="greater"))
		drawStar(homoVe3$p.value,5,6,2.2)#barely n.s.
		homoVnc <- with(interface.groups,wilcox.test(Homodimer,`SUMO-nc`,alternative="greater"))
		drawStar(homoVnc$p.value,2,6,2.65)

	},paste0(outdir,"interfaceBurial_",datalabel),9,6)

	# pdf(paste0(outdir,"accessibilityAndConservation.pdf"),7.5,4)
	html$figure(function(){

		# layout(rbind(1:2,3:4))
		op <- par(mfrow=c(2,3))
		xs <- boxplot(
			acc.groups,
			col=c("firebrick3","steelblue3","chartreuse3"),
			ylab="fitness score",xlab="Residue accessibility",
			ylim=c(-0.5,3.5)
		)
		coreSurf <- with(acc.groups,wilcox.test(core,surface,alternative="less"))
		if (coreSurf$p.value < 0.05) {
			y <- 3.4
			lines(c(1,1,3,3),c(y-.1,y,y,y-.1))
			text(2,y+.1,"*",cex=1.4)
		}
		ifSurf <- with(acc.groups,wilcox.test(interface,surface,alternative="less"))
		if (ifSurf$p.value < 0.05) {
			lines(c(2,2,3,3),c(3,3.1,3.1,3))
			text(2.5,3.2,"*",cex=1.4)
		}

		boxplot(
			cons.groups,
			col=c("gold1","orange","firebrick3"),
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

		boxplot(
			blosum.groups,
			names=c("< 0","≥ 0"),
			col=c("firebrick3","chartreuse3"),
			ylab="fitness score", xlab="BLOSUM62",
			ylim=c(-0.5,3.5)
		)
		blosum.test <- with(blosum.groups,wilcox.test(low,high,alternative="less"))
		if (blosum.test$p.value < 0.05) {
			y <- 3.4
			lines(c(1,1,2,2),c(y-.1,y,y,y-.1))
			text(1.5,y+.1,"*",cex=1.4)
		}

		boxplot(
			pp2.groups,
			col=c("chartreuse3","gold1","firebrick3"),
			names=c("benign","possibly\ndamaging","probably\ndamaging"),
			ylab="fitness score", xlab="Polyphen-2",
			ylim=c(-0.5,3.5)
		)
		hi.mid <- wilcox.test(pp2.groups[[1]],pp2.groups[[2]],alternative="greater")
		if (hi.mid$p.value < 0.05) {
			y <- 3.15
			lines(c(1,1,2,2),c(y-.1,y,y,y-.1))
			text(1.5,y+.1,"*",cex=1.4)
		}
		mid.lo <- wilcox.test(pp2.groups[[2]],pp2.groups[[3]],alternative="greater")
		if (mid.lo$p.value < 0.05) {
			y <- 3
			lines(c(2,2,3,3),c(y-.1,y,y,y-.1))
			text(2.5,y+.1,"*",cex=1.4)
		}

		boxplot(
			prov.groups,
			col=c("firebrick3","chartreuse3"),
			ylab="fitness score", xlab="PROVEAN",
			ylim=c(-0.5,3.5)
		)
		prov.test <- with(prov.groups,wilcox.test(Deleterious,Neutral),alternative="less")
		if (prov.test$p.value < 0.05) {
			y <- 3.4
			lines(c(1,1,2,2),c(y-.1,y,y,y-.1))
			text(1.5,y+.1,"*",cex=1.4)
		}

		boxplot(
			sift.groups,
			col=c("firebrick3","chartreuse3"),
			ylab="fitness score", xlab="SIFT",
			ylim=c(-0.5,3.5)
		)
		sift.test <- with(sift.groups,wilcox.test(Damaging,Tolerated),alternative="less")
		if (sift.test$p.value < 0.05) {
			y <- 3.4
			lines(c(1,1,2,2),c(y-.1,y,y,y-.1))
			text(1.5,y+.1,"*",cex=1.4)
		}

		par(op)

	},paste0(outdir,"accessibilityAndConservation_",datalabel),9,6)
	# invisible(dev.off())
}


data <- read.csv(paste0(outdir,"imputed_regularized_UBE2I_scores.csv"))
singles <- data[regexpr(",",data$mut) < 1,]
singles <- singles[regexpr("\\w\\d+\\w",singles$mut) > 0,]
singles <- singles[substr(singles$mut,1,1)!=substr(singles$mut,nchar(singles$mut),nchar(singles$mut)),]

calculateAndPlot(with(singles,data.frame(row.names=mut,mut=mut,score=joint.score)),"Regularized")

data <- read.csv(paste0(outdir,"compl_timeseries_results_byMut.csv"))
singles <- data[regexpr(",",data$mut) < 1,]
singles <- singles[regexpr("\\w\\d+\\w",singles$mut) > 0,]
singles <- singles[substr(singles$mut,1,1)!=substr(singles$mut,nchar(singles$mut),nchar(singles$mut)),]

calculateAndPlot(with(singles,data.frame(row.names=mut,mut=mut,score=score)),"BarSEQ")

data <- read.csv(paste0(outdir,"compl_tileSEQ_results_UBE2I_transformed.csv"))
singles <- data[regexpr(",|_",data$mut) < 1,]
singles <- singles[regexpr("\\w\\d+\\w",singles$mut) > 0,]
singles <- singles[substr(singles$mut,1,1)!=substr(singles$mut,nchar(singles$mut),nchar(singles$mut)),]

calculateAndPlot(with(singles,data.frame(row.names=mut,mut=mut,score=score.m)),"TileSEQ")

html$shutdown()

logger$info("Done.")
