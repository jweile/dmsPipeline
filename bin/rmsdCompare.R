
options(stringsAsFactors=FALSE)

source("lib/resultfile.R")
source("lib/liblogging.R")
source("lib/cliargs.R")

outdir <- getArg("outdir",default="workspace/test/")
# outdir <- "workspace/20170130-215857/"


#Initialize logger
logger <- new.logger(paste0(outdir,"rmsdCompare.log"))

#Set resultfile
html <- new.resultfile(paste0(outdir,"results.html"))
html$section("Error comparison")

geneNames <- c("UBE2I","SUMO1","TPK1","CALM1","NCS1")

doubleHist <- function(data1,data2,breaks=NULL,ymax=NULL,xlim=NULL,
		col=c("steelblue3","orange"),ylab1="data1",ylab2="data2",xlab="",
		main="") {

	hc <- hist(c(data1,data2),plot=FALSE)

	if (is.null(breaks)) {
		breaks <- hc$breaks
	}
	h1 <- hist(data1,plot=FALSE,breaks=breaks)
	h2 <- hist(data2,plot=FALSE,breaks=breaks)
	if (is.null(ymax)) {
		ymax <- max(h1$density,h2$density)
	}
	if (is.null(xlim)) {
		xlim <- range(breaks)
	}

	plot(0,type="n",xlim=xlim,ylim=c(-ymax,ymax),axes=FALSE,xlab=xlab,ylab="",main=main)
	axis(1)
	axis(2)
	mtext(ylab1,side=2,line=3,adj=1)
	mtext(ylab2,side=2,line=3,adj=0)

	rect(breaks[-length(breaks)],0,breaks[-1],h1$density,col=col[[1]])
	rect(breaks[-length(breaks)],0,breaks[-1],-h2$density,col=col[[2]])

}

infotable <- NULL

html$figure(function(){
	# layout(cbind(1:2,3:4,5:6,7:8,9:10))
	op <- par(mfrow=c(3,2))
	infotable <<- do.call(rbind,lapply(geneNames, function(geneName) {

		data <- read.csv(paste0(outdir,"imputed_regularized_",geneName,"_flipped_scores.csv"))
		#remove wt positions
		data <- data[with(data,substr(mut,1,1)!=substr(mut,nchar(mut),nchar(mut))),]

		screenSE <- data$screen.sd/sqrt(data$df)
		regSE <- with(data,joint.se[!is.na(screen.score)])

		breaks <- seq(0,4,0.01)
		doubleHist(screenSE,regSE,breaks=breaks,xlab="stderr",
			main=geneName,xlim=c(0,1),ylab1="experimental",ylab2="regularized")

		data.frame(
			gene=geneName,
			mutPossible=nrow(data),
			mutAchieved=sum(!is.na(data$screen.score)),
			mutAchievedPercent=100*sum(!is.na(data$screen.score))/nrow(data),
			rmsd=with(data,joint.se[which(is.na(screen.score))[[1]]]),
			maxSDExp=max(screenSE,na.rm=TRUE),
			maxSDReg=max(regSE,na.rm=TRUE)
		)

	}))
	par(op)

},paste0(outdir,"rmsdcompare"),7,15)

outfile <- paste0(outdir,"rmsdcompare.csv")
write.table(infotable,outfile,sep=",",quote=FALSE,row.names=FALSE)

html$shutdown()

logger$info("Done.")
