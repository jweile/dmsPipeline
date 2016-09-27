options(stringsAsFactors=FALSE)

library("EBImage")
#for some reason importing "methods" explicitly 
# is required for EBImage to work in script mode
library("methods")
source("lib/cliargs.R")
source("lib/resultfile.R")
source("lib/liblogging.R")


outdir <- getArg("outdir",default="workspace/test/")
imgdir <- getArg("imgdir",default="~/projects/spotting/complImg/final")

#Init logger
logger <- new.logger(paste0(outdir,"geneticInteractions.log"))

if (!file.exists(imgdir) || !(file.info(imgdir)$isdir)) {
	err <- paste(imgdir,"Is not a valid directory")
	logger$fatal(err)
	stop(err)
}

#Set resultfile
html <- new.resultfile(paste0(outdir,"results.html"))
html$section("Spotting assay result details")


drawFigure <- function(main.scores,filebase=imgdir) {

	#transforms image matrix to be compatible with R's image function
	my.image <- function(data,...) image(t(t(data)[nrow(t(data)):1,]),...)

	#target width and height of individual images
	w <- 112
	inh <- 32
	outh <- inh/2

	#load and prep individual image strips
	strips <- lapply(main.scores$id,function(id) {
		cat(id,"\t")
		if (id == "ctrl") {
			id <- "UBE2I-HYC-WT"
			ctrl <- TRUE
		} else ctrl <- FALSE
		#correct WT name
		if (id == "UBE2I-HYC-WT1") id <- "UBE2I-HYC-WT"
		#find correct image file for given clone
		filename <- list.files(filebase,pattern=id,full.names=TRUE)
		#if none exists, return empty image
		if (length(filename) == 0) return(matrix(1,w,outh))
		#if more than one exists, use the first
		if (length(filename) > 1) filename <- filename[[1]]
		
		#read the image
		img <- readImage(filename)
		#scale it to prespecified size
		img <- resize(img, w=w,h=inh)
		#transform to grayscale
		imgMatrix <- apply(imageData(img),1:2,mean)
		if (ctrl) {
			#extract 25C control
			 imgMatrix[,1:(outh-1)]
		} else {
			#cut away the 25C control
			imgMatrix[,outh:inh]
		}
	})
	cat("\n")
	#overwrite images with failed controls with empty matrices
	toDelete <- which(is.na(main.scores$spotting))
	for (i in toDelete) {
		strips[[i]] <- matrix(1,w,outh)
	}
	#join image matrices into one big matrix
	imgBlock <- t(do.call(cbind,strips))

	#define overlay colors
	redOver <- rgb(200,20,20,50,maxColorValue=255)
	greenOver <- rgb(20,200,20,50,maxColorValue=255)
	yellowOver <- rgb(200,200,20,50,maxColorValue=255)

	#set plot layout
	layout(rbind(1,2),heights=c(.6,.4))
	#set margins and label orientations for first subplot
	op <- par(mar=c(5,5,4,1)+.1,las=2)
	#draw barplot
	xs <- with(main.scores,barplot(
		score,
		ylim=c(
			min(c(-.5,score),na.rm=TRUE)-max(sd,na.rm=TRUE)/2,
			max(score,na.rm=TRUE)+max(sd,na.rm=TRUE)/2
		),
		ylab="screen score",
		names.arg=mut,
		xaxs="i",
		border="NA"
	))
	barwidth <- xs[2,1]-xs[1,1]
	with(main.scores,arrows(xs,score-sd/2,xs,score+sd/2,length=0.05,angle=90,code=3,col="gray40"))
	if (any(is.na(main.scores$score))) {
		is <- which(is.na(main.scores$score))
		rect(xs[is]-barwidth/2,-10,xs[is]+barwidth/2,10,density=10,border=NA,col="gray")
	}
	# if (main.scores$id[[1]]=="ctrl") {
		
	# }
	#add guide lines
	grid(NA,NULL)
	abline(h=c(0,1),col=c("firebrick3","chartreuse3"))
	#reset graphical parameters
	par(op)
	#set margins for second subplot
	op <- par(mar=c(3,5,0,1)+.1)
	#draw image matrix
	my.image(
		imgBlock,
		ylab="spotting assay\ndilutions",
		col=colorRampPalette(c("black","white"))(256),
		axes=FALSE,
		useRaster=TRUE
	)
	axis(2,at=seq(1/14,1-1/14,length.out=7),labels=7:1)
	#draw colored overlays
	rect(0,2/7,1,3/7,col=greenOver,border=NA)
	rect(0,0,1,2/7,col=yellowOver,border=NA)
	rect(0,6/7,1,1,col=redOver,border=NA)
	#reset graphical parameters
	par(op)

}


logger$info("Reading input")

imputed <- read.csv(paste0(outdir,"imputed_regularized_UBE2I_scores.csv"))
rownames(imputed) <- imputed$mut
barseq <- read.csv(paste0(outdir,"compl_timeseries_results_byClone.csv"))
rownames(barseq) <- barseq$id

#IMPORT SCORE MATRIX
spotting <- read.csv("input/spotting_clones_sanger.csv")

scores <- spotting
scores$score <- barseq[spotting$id,"score"]
scores[which(scores$id=="UBE2I-HYC-WT1"),"score"] <- 1
scores$sd <- barseq[spotting$id,"score.bsd"]
#sort by score
scores <- scores[order(scores$score,decreasing=TRUE),]

logger$info("Plotting spectrum clones")

#DRAW FIGURE FOR SPECTRUM CLONES
plotScores <- rbind(
	data.frame(category="control",id="ctrl",plate="ctrl",well=NA,mut="Yeast WT",score=NA,spotting=1,Sanger=TRUE,Sanger.result=NA,sd=NA),
	scores[scores$category == "control" & scores$Sanger & !is.na(scores$spotting),],
	# scores[!(scores$category %in% c("control","imputation")),]
	scores[scores$category == "spectrum" & scores$Sanger & !is.na(scores$spotting),]
)

# pdf("spottingSpectrum_clean.pdf",12,5)
html$subsection("Spectrum clones BarSEQ")
html$figure(function(){
	drawFigure(plotScores)
},paste0(outdir,"spottingSpectrum_barseq"),12,5)
# dev.off()


logger$info("Plotting fast clones")

#DRAW FIGURE FOR FAST CLONES
plotScores <- rbind(
	data.frame(category="control",id="ctrl",plate="ctrl",well=NA,mut="Yeast WT",score=NA,spotting=1,Sanger=TRUE,Sanger.result=NA,sd=NA),
	scores[scores$category == "control" & scores$Sanger & !is.na(scores$spotting),],
	scores[scores$category == "fast" & scores$Sanger & !is.na(scores$spotting),]
)

# pdf("spottingFast_clean.pdf",5,5)
html$subsection("Fast growing clones BarSEQ")
html$figure(function(){
	drawFigure(plotScores)
},paste0(outdir,"spottingFast_barseq"),5,5)
# dev.off()



scores <- spotting
scores$score <- imputed[spotting$mut,"joint.score"]
scores[which(scores$id=="UBE2I-HYC-WT1"),"score"] <- 1
scores$sd <- imputed[spotting$mut,"joint.sd"]
#sort by score
scores <- scores[order(scores$score,decreasing=TRUE),]

logger$info("Plotting spectrum clones")

#DRAW FIGURE FOR SPECTRUM CLONES
plotScores <- rbind(
	data.frame(category="control",id="ctrl",plate="ctrl",well=NA,mut="Yeast WT",score=NA,spotting=1,Sanger=TRUE,Sanger.result=NA,sd=NA),
	scores[scores$category == "control" & scores$Sanger & !is.na(scores$spotting),],
	# scores[!(scores$category %in% c("control","imputation")),]
	scores[scores$category == "spectrum" & scores$Sanger & !is.na(scores$spotting),]
)

# pdf("spottingSpectrum_clean.pdf",12,5)
html$subsection("Spectrum clones Regularized")
html$figure(function(){
	drawFigure(plotScores)
},paste0(outdir,"spottingSpectrum_regularized"),12,5)
# dev.off()


logger$info("Plotting fast clones")

#DRAW FIGURE FOR FAST CLONES
plotScores <- rbind(
	data.frame(category="control",id="ctrl",plate="ctrl",well=NA,mut="Yeast WT",score=NA,spotting=1,Sanger=TRUE,Sanger.result=NA,sd=NA),
	scores[scores$category == "control" & scores$Sanger & !is.na(scores$spotting),],
	scores[scores$category == "fast" & scores$Sanger & !is.na(scores$spotting),]
)

# pdf("spottingFast_clean.pdf",5,5)
html$subsection("Fast growing clones Regularized")
html$figure(function(){
	drawFigure(plotScores)
},paste0(outdir,"spottingFast_regularized"),5,5)
# dev.off()



logger$info("Plotting imputed clones")

#DRAW FIGURE FOR SPECTRUM CLONES
plotScores <- rbind(
	data.frame(category="control",id="ctrl",plate="ctrl",well=NA,mut="Yeast WT",score=NA,spotting=1,Sanger=TRUE,Sanger.result=NA,sd=NA),
	scores[scores$category == "control" & scores$Sanger & !is.na(scores$spotting),],
	# scores[!(scores$category %in% c("control","imputation")),]
	scores[scores$category == "imputation" & scores$Sanger & !is.na(scores$spotting),]
)

# pdf("spottingSpectrum_clean.pdf",12,5)
html$subsection("Imputation test set")
html$figure(function(){
	drawFigure(plotScores)
},paste0(outdir,"spottingImputation"),5,5)
# dev.off()


html$shutdown()

logger$info("Done!")
