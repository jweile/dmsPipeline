
###########################
# draw genophenogram plot #
###########################

# wt.aa = wildtype amino acid sequence
#
# List of mutations and their effects given in the following three vectors:
# pos = vector of positions
# mut.aa = vector of mutant AAs
# score = vector of scores
# error = vector of stderr values
# a = bezier transformation intensity (with -0.5 <= a <= 0.5)
# 
genophenogram <- function(wt.aa, pos, mut.aa, score, error=NULL, a=0) {

	bend <- function(x,a=0) {
		if (is.na(x)) return(NA)
		if (a == 0) return(x)
		if (abs(a) > 0.5) stop("parameter a must be between -0.5 and 0.5")
		if (x < 0) x <- 0
		if (x > 1) return(x)
		q <- (1 - 2*a) / (4*a)
		t <- ifelse(a >= 0,1,-1) * sqrt(x/(2*a) + q^2) - q
		(1+2*a)*t - 2*a*t^2
	}
	if (a != 0) {
		score <- sapply(score,bend,a=a)
	}

	layout(cbind(c(3,1),c(4,2)),widths=c(9.5,.5),heights=c(2,9))
	#Main plot
	###########
	aas <- c("A","V","L","I","M","F","Y","W","R","H","K","D","E","S","T","N","Q","G","C","P")
	op <- par(cex=.6,las=1,mar=c(5,4,0,0)+.1)
	plot(NA,type="n",
		xlim=c(0,length(wt.aa)+1),ylim=c(0,length(aas)+1),axes=FALSE,
		xlab="AA position",ylab="AA residue",main=""
	)
	axis(1,c(1,seq(5,length(wt.aa),5)))
	axis(2,at=1:20,labels=rev(aas))

	text(-1,c(16.5,8),c("hydrophobic","polar"),srt=90)
	text(-3,c(8.5,11),c("-","+"),srt=90)
	arrows(-2,c(3.6,12.6),-2,c(12.4,20.4),length=.02,angle=90,code=3)
	arrows(-4,c(7.6,9.6),-4,c(9.4,12.4),length=.02,angle=90,code=3)

	x <- pos
	y <- length(aas) - sapply(mut.aa,function(a)which(aas==a)) + 1
	colRamp <- colorRampPalette(c("royalblue3","white","firebrick3"))(11)
	colIdx <- sapply(score,function(s) {
		if (s <= 0) 1
		else if (s <= 0.2) 2
		else if (s <= 0.4) 3
		else if (s <= 0.6) 4
		else if (s <= .8) 5
		else if (s <= 1.2) 6
		else if (s <= 1.4) 7
		else if (s <= 1.6) 8
		else if (s <= 1.8) 9
		else if (s <= 2) 10
		else 11
	})
	cols <- colRamp[colIdx]

	#change wt positions to gold color
	cols[which(mut.aa==wt.aa[pos])] <- "lightgoldenrod1"

	rect(x-.5,y-.5,x+.5,y+.5,col=cols,border=NA)

	if (!is.null(error)) {
		e <- .5 * error / max(error)
		e[which(mut.aa==wt.aa[pos])] <- NA
		segments(x-e,y-e,x+e,y+e)
	}

	par(op)

	#####Legend
	###########
	op <- par(cex=.6,mar=c(5,0,0,4)+.1)
	plot(NA,type="n",xlim=c(-1,1),ylim=c(0,13),axes=FALSE,xlab="",ylab="")
	rect(0,0:11,1,1:12,col=c(colRamp[1:6],colRamp[6:11]),border=NA)
	rect(0,12,1,13,col="lightgoldenrod1",border=NA)
	if(!is.null(error)) {
		es <- seq(0,0.5,length.out=13)*13/21
		ys <- 0:12+0.5
		segments(-.5-es,ys-es,-.5+es,ys+es)
	}
	axis(4,at=c(.5,6,11.5,12.5),labels=c(0,1,2,"wt"))
	mtext("score",side=4,line=2,las=3,cex=0.7)
	if (!is.null(error)){
		me <- max(error)
		axis(2,at=c(.5,6.5,12.5),labels=signif(c(0,me/2,me),2))
		mtext("stderr",side=2,line=2,las=3,cex=0.7)
	}
	par(op)

	#Summary bars
	###########

	barvals <- do.call(rbind,tapply(colIdx,x,function(idxs) {
		table(factor(idxs,levels=1:11))
	}))/20
	barcums <- cbind(0,t(apply(barvals,1,function(x)sapply(1:11,function(i)sum(x[1:i])))))

	par(cex=.6,mar=c(0,4,1,0)+.1)
	plot(
		0,type="n",
		xlim=c(0,length(wt.aa)+1),
		ylim=c(0,1),
		axes=FALSE,xlab="",
		ylab="pos/neutral/neg"
	)
	n <- length(wt.aa)
	for (i in 1:11) {
		rect(1:n-.5,barcums[,i],1:n+.5,barcums[,i]+barvals[,i],col=colRamp[[i]],border=NA)	
	}
	axis(2)
	par(op)

}