
topoScatter <- function(x,y,resolution=20,thresh=5,
					topoCol=colorRampPalette(c("black","red","yellow"))(20), maxFreq=NULL,
					xlim=range(x,na.rm=TRUE),ylim=range(y,na.rm=TRUE),log="",...) {

	if (length(x) != length(y)) {
		stop("x and y must be of same length!")
	}

	#remove infinite values and na
	fin2 <- function(x) x[apply(x,1,function(.x) all(!is.na(.x) & is.finite(.x))),]
	xy <- cbind(x,y)
	xy <- fin2(xy)
	x <- xy[,1]
	y <- xy[,2]

	widenRange <- function(x) {
		a <- x[[1]]
		b <- x[[2]]
		.a <- a - (b-a)/10000
		.b <- b + (b-a)/10000
		if (a > 0 && .a < 0) .a <- a/10
		return(c(.a,.b))
	}
	xlim <- widenRange(xlim)
	ylim <- widenRange(ylim)

	xBins <- if (regexpr("x",log) > 0) {
		exp(log(10)*seq(log10(xlim[[1]]),log10(xlim[[2]]),length.out=resolution))
	} else {
		seq(xlim[[1]],xlim[[2]],length.out=resolution)
	}
	yBins <- if (regexpr("y",log) > 0) {
		exp(log(10)*seq(log10(ylim[[1]]),log10(ylim[[2]]),length.out=resolution))
	} else {
		seq(ylim[[1]],ylim[[2]],length.out=resolution)
	}
	bins <- matrix(0,ncol=length(xBins)-1,nrow=length(yBins)-1)

	for (i in 1:length(x)) {
		xbi <- max(which(xBins <= x[[i]]))
		ybi <- max(which(yBins <= y[[i]]))
		bins[[ybi,xbi]] <- bins[[ybi,xbi]]+1
	}

	restFilter <- sapply(1:length(x),function(i) {
		xbi <- max(which(xBins <= x[[i]]))
		ybi <- max(which(yBins <= y[[i]]))
		bins[[ybi,xbi]] <= thresh
	})
	xRest <- x[restFilter]
	yRest <- y[restFilter]

	maxFreq <- if (is.null(maxFreq)) max(bins) else {
		if (maxFreq < max(bins)) max(bins) else maxFreq
	}
	cat("maxFreq =",maxFreq,"\n")

	colMap <- apply(bins,c(1,2),function(v) {
		ci <- round((v/maxFreq) * (length(topoCol)-1) + 1)
		topoCol[[ci]]
	})

	plot(NULL,type="n",xlim=xlim,ylim=ylim,log=log,...)

	points(xRest,yRest,pch=20)
	for (xbi in 1:ncol(bins)) {
		for (ybi in 1:nrow(bins)) {
			if (bins[ybi,xbi] > thresh) {
				rect(xBins[[xbi]],yBins[[ybi]],xBins[[xbi+1]],yBins[[ybi+1]],col=colMap[ybi,xbi],border=NA)
			}
		}
	}

}

# x <- rnorm(10000,0,1)+10
# y <- rnorm(10000,0,1)+10
# topoScatter(x,y,resolution=60,xlab="foo",ylab="bar",thresh=2)
# topoScatter(x,y,log="xy",resolution=60,xlab="foo",ylab="bar")
