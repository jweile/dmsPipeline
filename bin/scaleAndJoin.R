####################################
# Transform TileSEQ data to BarSEQ #
# scale and join screen data       #
####################################

source("lib/resultfile.R")
source("lib/libyogitools.R")
source("lib/liblogging.R")
# source("lib/topoScatter.R")
source("lib/cliargs.R")
library("hash")

options(stringsAsFactors=FALSE)

#get output directory
outdir <- getArg("outdir",default="workspace/test/")

#Initialize logger
logger <- new.logger(paste0(outdir,"scaleAndJoin.log"))

#Set resultfile
html <- new.resultfile(paste0(outdir,"results.html"))
html$section("Re-scaling and Joining of datasets")

##############
# LOAD INPUT #
##############
logger$info("Reading input")

barUbe2iFile <- paste0(outdir,"compl_timeseries_results_byMut.csv")
tileUbe2iFile <- paste0(outdir,"compl_tileSEQ_results_UBE2I.csv")
tileSumo1File <- paste0(outdir,"compl_tileSEQ_results_SUMO1.csv")

ccbr <- read.csv(tileUbe2iFile)
rownames(ccbr) <- ccbr$mut
ltri <- read.csv(barUbe2iFile)
#reserve multimutants for later and exclude them for now
multi.muts <- ltri[regexpr(",",ltri$mut) > 0,]
ltri <- ltri[regexpr(",",ltri$mut) < 1,]
rownames(ltri) <- ltri$mut


#######################################
# Find transformation between screens #
#######################################
logger$info("Determining transformation function")

common.muts <- intersect(ccbr$mut,ltri$mut)
ltriVccbr <- data.frame(
	mut=common.muts,
	ltri.score=ltri[common.muts,"score"],
	ltri.sd=ltri[common.muts,"sd"],
	ccbr.score=ccbr[common.muts,"mean.lphi"],
	ccbr.sd=ccbr[common.muts,"bsd"]
)
splinetable <- with(ltriVccbr,data.frame(
	ltri=ltri.score,
	ccbr=ccbr.score,
	ccbr.exp=exp(ccbr.score)
))
splinetable <- fin2(splinetable)
z <- lm(splinetable[,1]~.,splinetable[,-1])
x <- seq(-2,1,.05)
y <- predict(z,data.frame(ccbr=x,ccbr.exp=exp(x)))

#################################
#Plot screens against each other#
#################################
logger$info(" -> Drawing plot")
# pdf(paste0(outdir,"scoreTransformationInput.pdf"))
html$subsection("Scaling Input")
html$figure(function(){
	plot(NA,type="n",
		xlim=c(-2,1),ylim=c(-.5,3),
		xlab="RegSEQ Screen Fitness Score",ylab="BarSEQ Screen Fitness Score"
	)
	with(ltriVccbr,arrows(ccbr.score-ccbr.sd/2,ltri.score,ccbr.score+ccbr.sd/2,ltri.score,length=.01,code=3,angle=90))
	with(ltriVccbr,arrows(ccbr.score,ltri.score-ltri.sd/2,ccbr.score,ltri.score+ltri.sd/2,length=.01,code=3,angle=90))
	abline(h=0:1,col=c("firebrick3","chartreuse3"))
	#add spline to plot
	lines(x,y,col="blue",lty="dashed",lwd=2)
},paste0(outdir,"scoreTransformationInput"))
# invisible(dev.off())


#################################
# Apply transformation to UBE2I #
#################################
logger$info("Applying transformation to test set")

#transform scores and stdevs
transf <- function(m,sd,z) {
	a <- coefficients(z)[["ccbr.exp"]]
	b <- coefficients(z)[["ccbr"]]
	ic <- coefficients(z)[["(Intercept)"]]
	c(
		m = a * exp(m) + b*m + ic,
		sd = sd * (a*exp(m) + b) 
	)
}

ccbr.trans <- t(apply(ltriVccbr[,-1],1,function(row) {
	transf(m=row[["ccbr.score"]],sd=row[["ccbr.sd"]],z=z)
}))
ltriVccbr <- cbind(ltriVccbr,ccbr.trans=ccbr.trans)


logger$info(" -> Drawing plot")

#Plot transformed CCBR vs LTRI with adjusted SD
# pdf(paste0(outdir,"scoreTransformationOutput.pdf"))
html$subsection("Scaling output")
html$figure(function(){
	plot(NA,type="n",
		xlim=c(-.5,3),ylim=c(-.5,3),
		xlab="Transformed RegSEQ Screen Fitness Score",ylab="BarSEQ Screen Fitness Score"
	)
	with(ltriVccbr,arrows(ccbr.trans.m-ccbr.trans.sd*5/2,ltri.score,ccbr.trans.m+ccbr.trans.sd*5/2,ltri.score,length=.01,code=3,angle=90))
	with(ltriVccbr,arrows(ccbr.trans.m,ltri.score-ltri.sd/2,ccbr.trans.m,ltri.score+ltri.sd/2,length=.01,code=3,angle=90))
	abline(h=0:1,v=0:1,col=c("firebrick3","chartreuse3"))
},paste0(outdir,"scoreTransformationOutput"))
# invisible(dev.off())


##########################################
# Apply transformation to all tileseq data #
##########################################
logger$info("Applying transformation to all TileSEQ")

#transform all remaining ccbr scores
ccbr.all.trans <- t(apply(ccbr[,-1],1,function(row) {
	transf(m=row[["mean.lphi"]],sd=row[["bsd"]],z=z)
}))
ccbr <- cbind(ccbr,score=ccbr.all.trans)

write.table(ccbr,paste0(outdir,"compl_tileSEQ_results_UBE2I_transformed.csv"),sep=",",row.names=FALSE)


#Load and transform SUMO data
sumo <- read.csv(tileSumo1File)
#transformation
sumo.trans <- t(apply(sumo[,-1],1,function(row) {
	transf(m=row[["mean.lphi"]],sd=row[["bsd"]],z=z)
}))
sumo <- cbind(sumo,score=sumo.trans)

#Write results for SUMO1
write.table(sumo,paste0(outdir,"compl_tileSEQ_results_SUMO1_transformed.csv"),sep=",",row.names=FALSE)


####################
# Joining datasets #
####################
logger$info("Joining TileSEQ and BarSEQ data")

#Build complete joint table
join.datapoints <- function(ms,sds) {
	ws <- (1/sds)/sum(1/sds)
	mj <- sum(ws*ms)
	vj <- sum(ws*(sds^2+ms^2)) -mj^2
	c(mj=mj,sj=sqrt(vj))
}

allmuts <- sort(union(ltri$mut,ccbr$mut))
joint.data <- as.df(lapply(allmuts, function(m) {
	#if the mutation is present in ltri set
	if (m %in% ltri$mut && !is.na(ltri[m,"score"]) && !is.infinite(ltri[m,"score"])) {
		#and also present in the ccbr set
		if (m %in% ccbr$mut && !is.na(ccbr[m,"score.m"]) && !is.infinite(ccbr[m,"score.m"])){#in both
			vs <- c(ltri[m,"score"],ccbr[m,"score.m"])
			sds <- c(ltri[m,"sd"],ccbr[m,"score.sd"]*5)
			joint <- join.datapoints(vs,sds)
			list(mut=m,score=joint[["mj"]],sd=joint[["sj"]])
		} else {#or only in ltri
			c(list(mut=m),ltri[m,c("score","sd"),drop=TRUE])
		}
	} else {
		#if it is only in the ccbr set
		if (m %in% ccbr$mut && !is.na(ccbr[m,"score.m"]) && !is.infinite(ccbr[m,"score.m"])){#only in ccbr
			list(mut=m,score=ccbr[m,"score.m"],sd=ccbr[m,"score.sd"]*5)
		} else {#if it's NA or infinite in both
			# cat("Excluding",m,"\n")
			list(mut=m,score=NA,sd=NA)
		}
	}
}))
rownames(joint.data) <- joint.data$mut
joint.data <- joint.data[!is.na(joint.data$score),]

#Re-add multimutants
joint.all <- rbind(joint.data,multi.muts[,-4])

#########################
# Write results to file #
#########################
logger$info("Joining TileSEQ and BarSEQ data")

write.table(joint.all,paste0(outdir,"compl_joint_results_UBE2I.csv"),sep=",",row.names=FALSE)

html$shutdown()

logger$info("Done.")
