options(stringsAsFactors=FALSE)

source("lib/resultfile.R")
source("lib/liblogging.R")
source("lib/cliargs.R")

outdir <- getArg("outdir",default="workspace/test/")
# outdir <- getArg("outdir",default="workspace/20170130-215857/")


#Initialize logger
logger <- new.logger(paste0(outdir,"invitro.log"))

#Set resultfile
html <- new.resultfile(paste0(outdir,"results.html"))
html$section("In vitro sumoylation comparison")


screen <- read.csv(paste0(outdir,"imputed_regularized_UBE2I_scores.csv"))
row.names(screen) <- screen$mut


invitro <- read.csv("res/bernier-villamor_sumoylation.csv")
colnames(invitro) <- sub("\\.",",",colnames(invitro))

invitro <- rbind(invitro,cbind(WT=c(1,0),sapply(colnames(invitro)[-1], function(mut) screen[mut,c("joint.score","joint.se")])))
missing <- which(is.na(invitro["joint.score",]))
invitro <- invitro[,-missing]

ord <- c(1,1+order(unlist(invitro["joint.score",-1]),decreasing=TRUE))
invitro2 <- as.matrix(invitro[,ord])


cp <- colorRampPalette(c("royalblue3","white","firebrick3"))(11)
col.fill <- function(vals) sapply(vals,function(s) {
	if (is.na(s)) "gray80" 
	else if (s <= 0) cp[[1]]
	else if (s <= 0.2) cp[[2]]
	else if (s <= 0.4) cp[[3]]
	else if (s <= 0.6) cp[[4]]
	else if (s <= .8) cp[[5]]
	else if (s <= 1.2) cp[[6]]
	else if (s <= 1.4) cp[[7]]
	else if (s <= 1.6) cp[[8]]
	else if (s <= 1.8) cp[[9]]
	else if (s <= 2) cp[[10]]
	else cp[[11]]
})

html$figure(function(){
	layout(rbind(1,2),heights=c(2,1))
	op <- par(las=2,mar=c(5,5,1,1)+.1)
	score <- unlist(invitro2["joint.score",])
	se <- unlist(invitro2["joint.se",])
	xs <- barplot(
		score,ylim=c(-0.5,2.5),
		col=col.fill(score),
		border="black",ylab="complementation"
	)
	arrows(
		xs[,1],score+se,xs[,1],score-se,
		code=3,angle=90,length=0.05,col="black"
	)
	abline(h=0:1,col="gray",lty="dotted")
	par(op)
	# abline(h=0:1,col=c("firebrick3","chartreuse3"))

	op <- par(mar=c(1,5,0,1)+.1,las=1)
	plot(NA,type="n",xlim=c(0,19),ylim=c(0,3),axes=FALSE,xlab="",ylab="")
	rect(0:18,2,1:19,3,col=col.fill(invitro2["RanGAP1",,drop=TRUE]),border=NA)
	rect(0:18,1,1:19,2,col=col.fill(invitro2["p53",,drop=TRUE]),border=NA)
	rect(0:18,0,1:19,1,col=col.fill(invitro2["IkBa",,drop=TRUE]),border=NA)
	rect(0,0,19,3,border="gray")
	axis(2,labels=c("IkBa","P53","RanGAP1"),at=1:3-0.5)
	par(op)
},paste0(outdir,"invitro"),7,4)



html$shutdown()

logger$info("Done.")
