######################################
# Evaluate all steps of the pipeline #
# against the spotting assay results #
######################################

options(stringsAsFactors=FALSE)

source("lib/resultfile.R")
source("lib/liblogging.R")
source("lib/cliargs.R")
source("lib/libroc.R")

outdir <- getArg("outdir",default="workspace/test/")

#Set resultfile
html <- new.resultfile(paste0(outdir,"results.html"))
html$section("Evaluate spotting assay")

#Init logger
logger <- new.logger(paste0(outdir,"evaluateSpotting.log"))

##############
# LOAD INPUT #
##############
logger$info("Reading input")

#
#Import data
#
spotting <- read.csv("input/spottingAssay_results.csv")
#remove duplicates, empties and incorrect clones
spotting <- spotting[!is.na(spotting$spotting),]
spotting <- spotting[!duplicated(spotting$id),]
spotting <- spotting[spotting$Sanger,]
rownames(spotting) <- spotting$id

spottingByMut <- with(spotting,tapply(spotting,mut,mean))

barcoded <- read.csv(paste0(outdir,"compl_timeseries_results_byClone.csv"))
rownames(barcoded) <- barcoded$id

regseq <- read.csv(paste0(outdir,"compl_tileSEQ_results_UBE2I.csv"))
rownames(regseq) <- regseq$mut

joint.imp <- read.csv(paste0(outdir,"imputed_regularized_UBE2I_scores.csv"))
joint.infVimp <- read.csv(paste0(outdir,"UBE2I_featureMatrix.csv"))[,c("mut","posAverageAvail","mmAverageAvail","giImputeAvail")]
if (!all(joint.imp$mut == joint.infVimp$mut)) {
	stop("Unable to join matrices!")
} else {
	joint.imp <- cbind(joint.imp,joint.infVimp[,2:4])
	rownames(joint.imp) <- joint.imp$mut
}
joint.only <- joint.imp[!is.na(joint.imp$screen.score),]




logger$info("Drawing correlation plots")

# pdf(paste0(outdir,"spottingEvalAll.pdf"),7,5)
html$subsection("Correlation with scores")
html$figure(function(){
	op <- par(mfrow=c(2,3))
	# BARCODED
	inner.ids <- intersect(spotting$id,barcoded$id)
	inner <- data.frame(
		spotting=spotting[inner.ids,"spotting"],
		screen=barcoded[inner.ids,"score"],
		screen.sd=barcoded[inner.ids,"score.bsd"],
		screen.se=barcoded[inner.ids,"score.bsd"]/sqrt(barcoded[inner.ids,"df"])
	)
	rownames(inner) <- inner.ids

	with(inner,{
		x <- jitter(spotting)
		plot(x,screen,pch=20,
			xlab="Spotting Assay Score",ylab="Screen Score",
			col="steelblue3",main="BarSeq Screen"
		)
		arrows(x,screen+screen.se,x,screen-screen.se,
			angle=90,length=0.01,code=3,col="steelblue3"
		)
		text(0.5,1.5,sprintf("SCC = %.02f",cor(spotting,screen,method="spearman")))
	})

	#REGSEQ
	inner.ids <- intersect(names(spottingByMut),regseq$mut)
	inner <- data.frame(
		spotting=spottingByMut[inner.ids],
		screen=regseq[inner.ids,"mean.lphi"],
		screen.sd=regseq[inner.ids,"bsd"],
		screen.se=regseq[inner.ids,"bsd"]/sqrt(regseq[inner.ids,"df"])
	)
	rownames(inner) <- inner.ids

	with(inner,{
		x <- jitter(spotting)
		plot(x,screen,pch=20,
			xlab="Spotting Assay Score",ylab="Screen Score",
			col="firebrick3",main="TileSeq Screen"
		)
		arrows(x,screen+screen.se,x,screen-screen.se,
			angle=90,length=0.01,code=3,col="firebrick3"
		)
		text(0.5,0,sprintf("SCC = %.02f",cor(spotting,screen,method="spearman")))
	})

	#JOINT SCORES
	inner.ids <- intersect(names(spottingByMut),joint.only$mut)
	inner <- data.frame(
		spotting=spottingByMut[inner.ids],
		screen=joint.only[inner.ids,"screen.score"],
		screen.sd=joint.only[inner.ids,"screen.sd"],
		screen.se=joint.only[inner.ids,"screen.sd"]/joint.only[inner.ids,"df"]
	)
	rownames(inner) <- inner.ids

	with(inner,{
		x <- jitter(spotting)
		plot(x,screen,pch=20,
			xlab="Spotting Assay Score",ylab="Joint Screen Score",
			col="purple",main="Joint Scores"
		)
		arrows(x,screen+screen.se,x,screen-screen.se,
			angle=90,length=0.01,code=3,col="purple"
		)
		text(0.5,1.5,sprintf("SCC = %.02f",cor(spotting,screen,method="spearman")))
	})

	#IMPUTED SCORES
	# imputed.only <- joint.imp[with(joint.imp,is.na(screen.score) & !(giImputeAvail | mmAverageAvail)),]
	impctrl <- with(spotting,mut[category=="imputation" & Sanger])
	imputed.only <- joint.imp[with(joint.imp,(mut %in% impctrl) & !(giImputeAvail | mmAverageAvail)),]
	inner.ids <- intersect(names(spottingByMut),imputed.only$mut)
	rmsd <- joint.imp[which(is.na(joint.imp$screen.score))[[1]],"joint.sd"]
	inner <- data.frame(
		spotting=spottingByMut[inner.ids],
		screen=imputed.only[inner.ids,"predicted.score"],
		screen.sd=rmsd,
		screen.se=rmsd
	)
	rownames(inner) <- inner.ids

	with(inner,{
		x <- jitter(spotting)
		plot(x,screen,pch=20,
			xlab="Spotting Assay Score",ylab="Imputed Score",
			col="orange",main="Imputed Scores"
		)
		arrows(x,screen+screen.se,x,screen-screen.se,
			angle=90,length=0.01,code=3,col="orange"
		)
		text(0.5,0.8,sprintf("SCC = %.02f",cor(spotting,screen,method="spearman")))
	})

	#ULTIMATE SCORES
	# inferred.only <- joint.imp[with(joint.imp,is.na(screen.score) & (posAverageAvail | mmAverageAvail)),]
	inner.ids <- intersect(names(spottingByMut),joint.imp$mut)
	inner <- data.frame(
		spotting=spottingByMut[inner.ids],
		screen=joint.imp[inner.ids,"joint.score"],
		screen.sd=joint.imp[inner.ids,"joint.sd"],
		screen.se=joint.imp[inner.ids,"joint.sd"]/joint.imp[inner.ids,"df"]
	)
	rownames(inner) <- inner.ids

	with(inner,{
		x <- jitter(spotting)
		plot(x,screen,pch=20,
			xlab="Spotting Assay Score",ylab="Regularized Score",
			col="chartreuse3",main="Regularized Scores"
		)
		arrows(x,screen+screen.se,x,screen-screen.se,
			angle=90,length=0.01,code=3,col="chartreuse3"
		)
		text(0.5,1.5,sprintf("SCC = %.02f",cor(spotting,screen,method="spearman")))
	})



	############
	# Completeness of each subset
	############
	completeness <- c(
		# BarSeq=length(unique(barcoded$mut[regexpr("^\\w\\d+\\w$",barcoded$mut) > 0]))/nrow(joint.imp),
		# TileSeq=length(unique(regseq$mut[regexpr("^\\w\\d+\\w$",regseq$mut) > 0 & is.finite(regseq$mean.lphi) & !is.na(regseq$mean.lphi)]))/nrow(joint.imp),
		BarSeq=length(unique(barcoded$mut[regexpr("^[A-Z]\\d+[A-Z]$",barcoded$mut) > 0]))/nrow(joint.imp),
		TileSeq=length(unique(regseq$mut[regexpr("^[A-Z]\\d+[A-Z]$",regseq$mut) > 0 & is.finite(regseq$mean.lphi) & !is.na(regseq$mean.lphi)]))/nrow(joint.imp),
		Joint=sum(!is.na(joint.imp$screen.score))/nrow(joint.imp),
		# Inferred=with(joint.imp,sum(is.na(screen.score) & (giImputeAvail | mmAverageAvail)))/sum(is.na(joint.imp$screen.score)),
		Imputed=1,
		Regularized=1
	)
	par(las=2,mar=c(6,4,4,1)+.1)
	barplot(completeness,ylab="completeness",main="Completeness")
	# par(op)

	par(op)
},paste0(outdir,"spottingEvalAll"),7,5)
# invisible(dev.off())


logger$info("Drawing ROC/PRC plots")


###################
#ROC & PRC CURVES #
###################

html$subsection("ROC and PRC curves")
html$figure(function(){
	op <- par(mfrow=c(4,2))
	# BARCODED
	inner.ids <- intersect(spotting$id,barcoded$id)
	inner <- data.frame(
		spotting=spotting[inner.ids,"spotting"],
		screen=barcoded[inner.ids,"score"],
		screen.sd=barcoded[inner.ids,"score.bsd"]
	)
	rownames(inner) <- inner.ids

	rd <- roc.data(truth=(inner$spotting >= 0.75),scores=inner$screen)
	draw.roc(rd,main="BarSeq")
	text(50,50,sprintf("AUROC = %.02f",auroc(rd)))
	draw.prc(rd)
	text(50,50,sprintf("AUPRC = %.02f",auprc(rd)))


	#REGSEQ
	inner.ids <- intersect(names(spottingByMut),regseq$mut)
	inner <- data.frame(
		spotting=spottingByMut[inner.ids],
		screen=regseq[inner.ids,"mean.lphi"],
		screen.sd=regseq[inner.ids,"bsd"]*5
	)
	rownames(inner) <- inner.ids

	rd <- roc.data(truth=(inner$spotting >= 0.75),scores=inner$screen)
	draw.roc(rd,main="TileSeq")
	text(50,50,sprintf("AUROC = %.02f",auroc(rd)))
	draw.prc(rd)
	text(50,50,sprintf("AUPRC = %.02f",auprc(rd)))

	#JOINT SCORES
	inner.ids <- intersect(names(spottingByMut),joint.only$mut)
	inner <- data.frame(
		spotting=spottingByMut[inner.ids],
		screen=joint.only[inner.ids,"screen.score"],
		screen.sd=joint.only[inner.ids,"screen.sd"]
	)
	rownames(inner) <- inner.ids

	rd <- roc.data(truth=(inner$spotting >= 0.75),scores=inner$screen)
	draw.roc(rd,main="Joint screens")
	text(50,50,sprintf("AUROC = %.02f",auroc(rd)))
	draw.prc(rd)
	text(50,50,sprintf("AUPRC = %.02f",auprc(rd)))

	#INFERRED SCORES
	# inferred.only <- joint.imp[with(joint.imp,is.na(screen.score) & (giImputeAvail | mmAverageAvail)),]
	# inner.ids <- intersect(names(spottingByMut),inferred.only$mut)
	# inner <- data.frame(
	# 	spotting=spottingByMut[inner.ids],
	# 	screen=inferred.only[inner.ids,"predicted.score"],
	# 	screen.sd=rep(0.4285,length(inner.ids))
	# )
	# rownames(inner) <- inner.ids

	# rd <- roc.data(truth=(inner$spotting >= 0.75),scores=inner$screen)
	# draw.roc(rd,main="Inferred")
	# text(50,50,sprintf("AUROC = %.02f",auroc(rd)))
	# draw.prc(rd)
	# text(50,50,sprintf("AUPRC = %.02f",auprc(rd)))

	#ULTIMATE SCORES
	# inferred.only <- joint.imp[with(joint.imp,is.na(screen.score) & (posAverageAvail | mmAverageAvail)),]
	inner.ids <- intersect(names(spottingByMut),joint.imp$mut)
	inner <- data.frame(
		spotting=spottingByMut[inner.ids],
		screen=joint.imp[inner.ids,"joint.score"],
		screen.sd=joint.imp[inner.ids,"joint.sd"]
	)
	rownames(inner) <- inner.ids

	rd <- roc.data(truth=(inner$spotting >= 0.75),scores=inner$screen)
	draw.roc(rd,main="Regularized")
	text(50,50,sprintf("AUROC = %.02f",auroc(rd)))
	draw.prc(rd)
	text(50,50,sprintf("AUPRC = %.02f",auprc(rd)))

	par(op)
},paste0(outdir,"spottingROC"),4,10)

html$shutdown()

logger$info("Done.")
