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

ccbr$se <- with(ccbr,bsd/sqrt(df))


#######################################
# Find transformation between screens #
#######################################
logger$info("Determining transformation function")

common.muts <- intersect(ccbr$mut,ltri$mut)
ltriVccbr <- data.frame(
	mut=common.muts,
	ltri.score=ltri[common.muts,"score"],
	ltri.sd=ltri[common.muts,"sd"],
	ltri.se=ltri[common.muts,"se"],
	ltri.df=ltri[common.muts,"df"],
	ccbr.score=ccbr[common.muts,"mean.lphi"],
	ccbr.sd=ccbr[common.muts,"bsd"],
	ccbr.se=ccbr[common.muts,"se"],
	ccbr.df=ccbr[common.muts,"df"]
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
		xlim=c(-2,1),ylim=c(-.5,2.5),
		xlab="RegSEQ Screen Fitness Score",ylab="BarSEQ Screen Fitness Score"
	)
	with(ltriVccbr,arrows(ccbr.score-ccbr.se,ltri.score,ccbr.score+ccbr.se,ltri.score,length=.01,code=3,angle=90))
	with(ltriVccbr,arrows(ccbr.score,ltri.score-ltri.se,ccbr.score,ltri.score+ltri.se,length=.01,code=3,angle=90))
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
ltriVccbr$ccbr.trans.se <- with(ltriVccbr,ccbr.trans.sd/sqrt(ccbr.df))


logger$info(" -> Drawing plot")

#Plot transformed CCBR vs LTRI with adjusted SD
# pdf(paste0(outdir,"scoreTransformationOutput.pdf"))
html$subsection("Scaling output")
html$figure(function(){
	plot(NA,type="n",
		xlim=c(-.5,2.5),ylim=c(-.5,2.5),
		xlab="Transformed RegSEQ Screen Fitness Score",ylab="BarSEQ Screen Fitness Score"
	)
	with(ltriVccbr,arrows(ccbr.trans.m-ccbr.trans.se,ltri.score,ccbr.trans.m+ccbr.trans.se,ltri.score,length=.01,code=3,angle=90))
	with(ltriVccbr,arrows(ccbr.trans.m,ltri.score-ltri.se,ccbr.trans.m,ltri.score+ltri.se,length=.01,code=3,angle=90))
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
ccbr$score.se <- with(ccbr,score.sd/sqrt(df))

ub.outfile <- paste0(outdir,"compl_tileSEQ_results_UBE2I_transformed.csv")
write.table(ccbr,ub.outfile,sep=",",row.names=FALSE)


#Load and transform SUMO data
sumo <- read.csv(tileSumo1File)
#transformation
sumo.trans <- t(apply(sumo[,-1],1,function(row) {
	transf(m=row[["mean.lphi"]],sd=row[["bsd"]],z=z)
}))
sumo <- cbind(sumo,score=sumo.trans)
sumo$score.se <- with(sumo,score.sd/sqrt(df))


#Write results for SUMO1
sumo.outfile <- paste0(outdir,"compl_tileSEQ_results_SUMO1_transformed.csv")
write.table(sumo,sumo.outfile,sep=",",row.names=FALSE)

html$subsection("Transformed data")
html$link.data(ub.outfile)
html$link.data(sumo.outfile)


####################
# Joining datasets #
####################
logger$info("Joining TileSEQ and BarSEQ data")

#Build complete joint table
join.datapoints <- function(ms,sds,ses) {
	ws <- (1/ses)/sum(1/ses)
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
			sds <- c(ltri[m,"sd"],ccbr[m,"score.sd"])
			ses <- c(ltri[m,"se"],4*ccbr[m,"score.se"])
			dfs <- ltri[m,"df"]+ccbr[m,"df"]
			joint <- join.datapoints(vs,sds,ses)
			list(mut=m,score=joint[["mj"]],sd=joint[["sj"]],df=dfs)
		} else {#or only in ltri
			c(list(mut=m),ltri[m,c("score","sd","df"),drop=TRUE])
		}
	} else {
		#if it is only in the ccbr set
		if (m %in% ccbr$mut && !is.na(ccbr[m,"score.m"]) && !is.infinite(ccbr[m,"score.m"])){#only in ccbr
			list(mut=m,score=ccbr[m,"score.m"],sd=ccbr[m,"score.sd"],df=ccbr[m,"df"])
		} else {#if it's NA or infinite in both
			# cat("Excluding",m,"\n")
			list(mut=m,score=NA,sd=NA,df=NA)
		}
	}
}))
rownames(joint.data) <- joint.data$mut
joint.data <- joint.data[!is.na(joint.data$score),]

#Re-add multimutants
joint.all <- rbind(joint.data,multi.muts[,-4])

joint.all$se <- with(joint.all,sd/sqrt(df))

#########################
# Write results to file #
#########################
logger$info("Joining TileSEQ and BarSEQ data")

outfile <- paste0(outdir,"compl_joint_results_UBE2I.csv")
write.table(joint.all,outfile,sep=",",row.names=FALSE)

html$subsection("Joined data")
html$link.dat(outfile)

html$shutdown()

logger$info("Done.")
