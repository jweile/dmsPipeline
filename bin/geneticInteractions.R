options(stringsAsFactors=FALSE)


source("lib/resultfile.R")
source("lib/libyogitools.R")
source("lib/liblogging.R")
source("lib/cliargs.R")
source("lib/topoScatter.R")

library("hash")


#Get output directory
outdir <- getArg("outdir",default="workspace/test/")
infile <- getArg("infile",default="workspace/test/compl_joint_results_UBE2I.csv")

#Set resultfile
html <- new.resultfile(paste0(outdir,"results.html"))
html$section("Intragenic epistasis")

#Init logger
logger <- new.logger(paste0(outdir,"geneticInteractions.log"))

##############
# LOAD INPUT #
##############
logger$info("Reading input")

data <- read.csv(infile)
muts <- strsplit(data$mut,",")
data$single <- sapply(muts,length) < 2
data$double <- sapply(muts,length) == 2


#########################
# BUILD CANDIDATE TABLE #
#########################

logger$info("Indexing single mutants")

single.idx <- hash()
for (i in 1:nrow(data)) {
	if (data$single[[i]]) {
		m <- data$mut[[i]]
		single.idx[[m]] <- i
	}
}

logger$info("Building double mutant table")

dm.data <- with(data,to.df(do.call(rbind,lapply(which(double),function(i) {
	ms <- muts[[i]]
	i1 <- single.idx[[ms[[1]]]]
	i2 <- single.idx[[ms[[2]]]]
	if (is.null(i1) || is.null(i2)) {
		return(rep(NULL,7))
	}
	list(
		dm=mut[[i]],dm.fitness=score[[i]],dm.fitness.sd=sd[[i]],
		sm1.fitness=score[[i1]], sm1.fitness.sd=sd[[i1]],
		sm2.fitness=score[[i2]], sm2.fitness.sd=sd[[i2]]
	)
}))))

logger$info("Correcting negative fitness values")

zero.bound <- function(xs) sapply(xs,function(x) ifelse(x<0,0,x))
dm.data2 <- within(dm.data,{
	dm.fitness <- zero.bound(dm.fitness)
	sm1.fitness <- zero.bound(sm1.fitness)
	sm2.fitness <- zero.bound(sm2.fitness)
})

#####################################
# CALCULATE EXPECTED FITNESS AND SD #
#####################################

logger$info("Calculating expected fitness and stdev")

prod.sd <- function(e1,sd1,e2,sd2) sqrt(sd1^2 * sd2^2 + sd1^2 * e2^2 + sd2^2 * e1^2)

dm.data2$edmf <- with(dm.data2, sm1.fitness*sm2.fitness)
dm.data2$edmf.sd <- with(dm.data2,mapply(prod.sd,
	e1=sm1.fitness,sd1=sm1.fitness.sd,e2=sm2.fitness,sd2=sm2.fitness.sd
))
dm.data2$epsilon <- with(dm.data2,dm.fitness-edmf)

#########################
# SIGNIFICANCE TESTING  #
#########################

logger$info("Performing t-tests and FDR correction")

tstat <- function (n1,n2,m1,m2,sd1,sd2) {
  tt <- -(m1 - m2)/sqrt((((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2)/(n1 + n2 - 2)) * ((n1 + n2)/(n1 * n2)))
  dft <- n1 + n2 - 2
  # t^2 follows F distribution
  p <- 1 - pf(tt^2, 1, dft)
  c(t=tt,p=p)
}

t.result <- with(dm.data2,to.df(do.call(rbind,mapply(tstat,
	n1=6,n2=6,m1=dm.fitness,m2=edmf,sd1=dm.fitness.sd,sd2=edmf.sd,
	SIMPLIFY=FALSE
))))
t.result$q <- p.adjust(t.result$p,method="fdr")

dm.data2 <- cbind(dm.data2,t.result)


logger$info("Filtering for significance and sorting by effect size")

gis <- dm.data2[dm.data2$q < .05,]
gis <- gis[order(abs(gis$epsilon),decreasing=TRUE),]


#################
# VOLCANO PLOT  #
#################

#Volcano plot
logger$info("Drawing volcano plot")

# pdf(paste0(outdir,"epistasis_volcano.pdf"),5,5)
html$figure(function(){
	layout(rbind(1,2),heights=c(2,5))
	op <- par(mar=c(0,4,4,1)+.1,las=1)
	with(dm.data2,hist(
		epsilon,breaks=20,
		prob=TRUE,
		xlab="",axes=FALSE,main="UBE2I Intragenic Epistasis",
		col="gray",border=NA
	))
	axis(2)
	par(mar=c(5,4,0,1)+.1)
	with(dm.data2,plot(epsilon,-log10(q),
		pch=16,	col=rgb(79,148,205,80,maxColorValue=255),
		# resolution=60, thresh=3,
		xlab=expression(epsilon),ylab=expression(-log[10](q))
	))
	abline(h=-log10(.05),col="red",lty="dashed")
	text(-1.7,-log10(.05),"q = 0.05",col="red",pos=3)
	par(op)
}, paste0(outdir,"epistasis_volcano"),5,5)
# invisible(dev.off())


################
# MATRIX PLOT  #
################

gimat <- matrix(0,nrow=160,ncol=160)
for (i in 1:nrow(gis)) {
	muts <- strsplit(gis$dm[[i]],",")[[1]]
	pos <- as.integer(substr(muts,2,nchar(muts)-1))
	gimat[pos[[1]],pos[[2]]] <- gis$epsilon[[i]]
	gimat[pos[[2]],pos[[1]]] <- gis$epsilon[[i]]
}

html$figure(function(){
	layout(cbind(1,2),widths=c(8,2))
	op <- par(mar=c(5,4,4,0)+.1)
	colramp <- colorRampPalette(c("blue","white","red"))(15)
	image(
		gimat,
		zlim=c(-2.5,2.5),
		col=colramp,
		axes=FALSE,xlab="AA pos",ylab="AA pos"
	)
	#2-35, 36-71, 72-105, 106-142, 143-158
	starts <- c(2,36,72,106,143)/160
	stops <- c(35,71,105,142,158)/160
	rect(starts,starts,stops,stops,border="gray",lty="dashed")
	ticks <- c(1,seq(10,160,10))
	axis(1,at=ticks/160,labels=ticks)
	axis(2,at=ticks/160,labels=ticks)
	op <- par(mar=c(5,1,4,4)+.1)
	plot(0,type="n",axes=FALSE,xlim=c(0,1),ylim=c(0,15),xlab="",ylab="")
	rect(0,0:14,1,1:15,col=colramp,border=NA)
	axis(4,at=c(0.5,7.5,14.5),labels=c(-2.5,0,2.5))
	mtext(expression(epsilon),side=4,line=3)
	par(op)
}, paste0(outdir,"epistasis_matrix"),7,7)


###########
# OUTPUT  #
###########

logger$info("Writing results to file")
outfile <- paste0(outdir,"genetic_interactions.csv")
write.table(gis,outfile,sep=",",row.names=FALSE)
html$subsection("Output")
html$link.data(outfile)

html$shutdown()

logger$info("Done.")
