# features <- read.csv("features.csv",stringsAsFactors=FALSE)
# features <- read.csv("burialsecstruc.csv",stringsAsFactors=FALSE)

# plot(NA,xlim=c(0,10),ylim=c(0,2))
# draw.helix(x=c(0,5),.5)
# draw.sheet(x=c(6,10),.5)

new.trackdrawer <- function(l) {

	applyCol <- function(vals, cols) sapply(vals,function(x) if(is.na(x)) "gray" else cols[round(x*10+1)])

	.tracks <- list()
	.labels <- character()
	.types <- character()

	add.track <- function(values, label, col, minVal=0,maxVal=max(values,na.rm=TRUE)) {
		colRamp <- colorRampPalette(c("white",col))(11)
		normvals <- (values - minVal)/(maxVal-minVal)
		colVals <- applyCol(normvals,colRamp)
		.tracks[[length(.tracks)+1]] <<- colVals
		.labels[[length(.labels)+1]] <<- label
		.types[[length(.types)+1]] <<- "heat"
	}

	draw.helix <- function(xs,y) {
		loop <- cbind(
			x=c(0,.1,.2,.3,.4,.5,.6,.7,.75,.7,.6,.5,.4,.3,.25,.3,.4,.5,.6,.7,.8,.9,1),
			y=c(.1,.11,.13,.15,.17,.19,.24,.32,.5,.75,.87,.9,.87,.75,.5,.32,.24,.19,.17,.15,.13,.11,.1)
		)
		draw.loop <- function(x,y) {
			offset <- cbind(rep(x,nrow(loop)),rep(y,nrow(loop)))
			lines(loop+offset,lwd=1)
		}
		for (x in min(xs):max(xs)) draw.loop(x,y)
	}

	draw.sheet <- function(xs,y) {
		start <- min(xs)
		end <- max(xs)
		tipend <- end-.8
		polygon(
			x=c(start,tipend,tipend,end,tipend,tipend,start),
			y=y+c(.2,.2,0,.5,1,.8,.8),
			lwd=1
		)
	}

	draw.dis <- function(xs,y) {
		start <- min(xs)
		end <- max(xs)
		arrows(start,y+.5,end,y+.5,code=0,lty="dashed")
	}

	add.ss.track <- function(values) {
		start <- -1
		track <- list()
		for (i in 1:length(values)) {
			if (!is.na(values[[i]]) && values[[i]] == "AlphaHelix") {
				#if this is the start
				if (i==1 || is.na(values[[i-1]]) || values[[i-1]] != "AlphaHelix") {
					#if previous was strand, then end strand
					if (i > 1 && !is.na(values[[i-1]]) && values[[i-1]] == "Strand") {
						track[[length(track)+1]] <- list(f=draw.sheet,xs=c(start,i-1))
					#if previous was disorder, then end disorder
					} else if (i > 1 && !is.na(values[[i-1]]) && values[[i-1]] == "Disorder") {
						track[[length(track)+1]] <- list(f=draw.dis,xs=c(start,i-1))
					}
					start <- i
				}
			} else if (!is.na(values[[i]]) && values[[i]] == "Strand") {
				#if this is the start
				if (i==1 || is.na(values[[i-1]]) || values[[i-1]] != "Strand") {
					#if previous was helix, then end strand
					if (i > 1 && !is.na(values[[i-1]]) && !is.na(values[[i-1]]) && values[[i-1]] == "AlphaHelix") {
						track[[length(track)+1]] <- list(f=draw.helix,xs=c(start,i-1))
					#if previous was disorder, then end strand
					} else if (i > 1 && !is.na(values[[i-1]]) && values[[i-1]] == "Disorder") {
						track[[length(track)+1]] <- list(f=draw.dis,xs=c(start,i-1))
					}
					start <- i
				}
			} else if (!is.na(values[[i]]) && values[[i]] == "Disorder") {
				#if this is the start
				if (i==1 || is.na(values[[i-1]]) || values[[i-1]] != "Disorder") {
					#if previous was helix, then end strand
					if (i > 1 && !is.na(values[[i-1]]) && values[[i-1]] == "AlphaHelix") {
						track[[length(track)+1]] <- list(f=draw.helix,xs=c(start,i-1))
					#if previous was strand, then end strand
					} else if (i > 1 && !is.na(values[[i-1]]) && values[[i-1]] == "Strand") {
						track[[length(track)+1]] <- list(f=draw.sheet,xs=c(start,i-1))
					}
					start <- i
				}
			} else {
				if (i > 1 && !is.na(values[[i-1]]) && values[[i-1]] == "AlphaHelix") {
					track[[length(track)+1]] <- list(f=draw.helix,xs=c(start,i-1))
				} else if (i > 1 && !is.na(values[[i-1]]) && values[[i-1]] == "Strand") {
					track[[length(track)+1]] <- list(f=draw.sheet,xs=c(start,i-1))
				} else if (i > 1 && !is.na(values[[i-1]]) && values[[i-1]] == "Disorder") {
					track[[length(track)+1]] <- list(f=draw.dis,xs=c(start,i-1))
				}
			}
		}
		.tracks[[length(.tracks)+1]] <<- track
		.labels[[length(.labels)+1]] <<- ""
		.types[[length(.types)+1]] <<- "ss"
	}

	draw <- function() {
		if (length(.tracks) < 1) stop("no tracks!")
		h <- length(.tracks)
		op <- par(las=1,mar=c(5,7,4,2)+.1)
		plot(NA,type="n",xlim=c(0,l),ylim=c(0,h),
			xlab="AA position",axes=FALSE,ylab=""
		)
		axis(1,at=c(1,seq(5,l,5))-.5,labels=c(1,seq(5,l,5)))
		axis(2,at=h:1-.5,labels=.labels)

		for (i in 1:h) {
			if(.types[[i]] == "ss") {
				for (elem in .tracks[[i]]) {
					elem$f(elem$xs,y=h-i)
				}
			} else if (.types[[i]] == "heat") {
				pos <- 1:l
				rect(pos-1,h-i,pos,h-i+1,border=NA,col=.tracks[[i]])
			}
		}
		par(op)
	}

	debug <- function() {
		cat("Tracks:\n")
		print(.tracks)
		cat("Labels:\n")
		print(.labels)
		cat("Types:\n")
		print(.types)
	}

	list(add.track=add.track,add.ss.track=add.ss.track,draw=draw,debug=debug)
}

# trackdrawer <- new.trackdrawer(l=nrow(features))
# trackdrawer$add.ss.track(features$ssName)
# trackdrawer$add.track(features$ratio.accessible,"Accessibility","steelblue3")
# trackdrawer$add.track(features$senp1,"SENP1","orange",maxVal=1)
# trackdrawer$add.track(features$senp2,"SENP2","orange",maxVal=1)
# trackdrawer$add.track(features$e1,"UBA2 (E1)","orange",maxVal=1)
# trackdrawer$add.track(features$ube2i_cov,"UBE2I (E2)","orange",maxVal=1)
# trackdrawer$add.track(features$ube2i_noncov,"UBE2I (E2) NC","orange",maxVal=1)
# trackdrawer$add.track(features$ranbp,"RanBP2 (E3)","orange",maxVal=1)
# trackdrawer$add.track(features$daxx,"DAXX (SIM)","orange",maxVal=1)
# trackdrawer$add.track(features$pias,"PIAS (SIM)","orange",maxVal=1)
# trackdrawer$add.track(features$pml,"PML (SIM)","orange",maxVal=1)
# trackdrawer$add.track(features$tdg,"TDG (SIM)","orange",maxVal=1)
# trackdrawer$draw()

# pdf("burial_tracks.pdf",10,4)
# trackdrawer$draw()
# dev.off()