####################################
# Process Raw mutation counts from #
# TileSEQ data                     #
####################################

source("lib/resultfile.R")
source("lib/libyogitools.R")
source("lib/liblogging.R")
source("lib/topoScatter.R")
source("lib/cliargs.R")
library("hash")

options(stringsAsFactors=FALSE)

#get output directory
outdir <- getArg("outdir",default="workspace/test/")
#get input file
infile <- getArg("infile",default="input/raw_counts_UBE2I_tileseq.tsv")
#get gene name
geneName <- getArg("geneName",default="UBE2I")

#Initialize logger
logger <- new.logger(paste0(outdir,"tileseq-",geneName,".log"))

#Set resultfile
html <- new.resultfile(paste0(outdir,"results.html"))
html$section(paste(geneName,"Complementation TileSEQ Data Analysis"))

##############
# LOAD INPUT #
##############
logger$info("Reading input")

ccbr <- read.delim(infile)
ccbr.mut <- with(ccbr,paste(wt_aa,pos,mut_aa,sep=""))

#extract condition names
conditions <- colnames(ccbr)[-c(1:6)]

##############
# PROCESS #
##############
logger$info("Processing data")

#collapse codons into unique AA changes
ccbr2 <- as.df(tapply(1:nrow(ccbr),ccbr.mut,function(is) {
	mut <- unique(ccbr.mut[is])
	c(list(mut=mut),colSums(ccbr[is,conditions]))
},simplify=FALSE))

#remove rows with zeroes
ccbr2 <- ccbr2[apply(ccbr2[,conditions[1:2]],1,function(vs) !any(vs==0)),]
#remove rows with mean(controls) + 3SD > min(mean(select),mean(nonselect))
flagged <- apply(ccbr2[,conditions],1,function(vs) {
	ctr <- mean(vs[c("control1","control2")])
	sdctr <- sd(vs[c("control1","control2")])
	sel <- mean(vs[c("select1","select2")])
	nsel <- mean(vs[c("nonselect1","nonselect2")])
	ctr + 3*sdctr > min(sel,nsel)
})
ccbr2 <- ccbr2[!flagged,]

#Calculate fitness phi replicates and error metrics
ccbr3 <- with(ccbr2,data.frame(
	mut=mut,
	phi1=(select1 - control1) / (nonselect1 - control1),
	phi2=(select2 - control2) / (nonselect2 - control2),
	minBC=apply(cbind(nonselect1,nonselect2),1,min)
))
ccbr3$sd <- apply(ccbr3[,c("phi1","phi2")],1,sd)
ccbr3$meanphi <-  apply(ccbr3[,c("phi1","phi2")],1,mean)



logger$info("Plotting replicate correlation")

#Correlation Dot-plot
# pdf(paste0(outdir,"tileSEQ_",geneName,"_phi_replicates.pdf"),5,5)
html$subsection("Replicate correlation")
html$figure(function(){
	with(ccbr3[ccbr3$minBC > 500,],topoScatter(phi1,phi2,
		xlab=expression(phi[1]),ylab=expression(phi[2]),
		resolution=50
	))
},paste0(outdir,"tileSEQ_",geneName,"_phi_replicates"),5,5)
# invisible(dev.off())


logger$info("Plotting Error Regularization W/O log transformation")

# pdf(paste0(outdir,"tileSEQ_",geneName,"_SDvPHI.pdf"),10,5)
html$subsection("Modeling Error Prior without log transformation")
html$figure(function(){
	op <- par(mfrow=c(1,2))
	#Plot minBC vs SD
	with(ccbr3,topoScatter(minBC,sd+0.00001,log="xy",maxFreq=35,thresh=3,
		resolution=40, xlab="Read depth (Millions)", ylab=expression(sigma)
	))
	#Plot score vs SD
	with(ccbr3[ccbr3$minBC > 500,],topoScatter(meanphi,sd,log="xy",pch=20,resolution=40,
		xlab=expression(E(phi)),ylab=expression(sigma),maxFreq=35,thresh=3
	))
	par(op)
},paste0(outdir,"tileSEQ_",geneName,"_SDvPHI"),10,5)
# invisible(dev.off())


logger$info("Log transforming data")

ccbr.log <- with(ccbr3,data.frame(
	mut=mut,
	lphi1=log10(phi1),
	lphi2=log10(phi2),
	sd=apply(log10(cbind(phi1,phi2)),1,sd)
))
ccbr.log$mean.lphi <- apply(ccbr.log[,c("lphi1","lphi2")],1,mean)


logger$info("Plotting Error Regularization WITH log transformation")

# pdf(paste0(outdir,"tileSEQ_",geneName,"_errorRegularizationInput.pdf"),10,5)
html$subsection("Modeling Error Prior with log transformation")
html$figure(function(){
	op <- par(mfrow=c(1,2))
	#Plot minBC vs SD
	with(ccbr.log,topoScatter(ccbr3$minBC,sd+0.00001,log="xy",maxFreq=35,thresh=3,
		resolution=40, xlab="Read depth (Millions)", ylab=expression(sigma)
	))
	#Plot score vs SD
	with(ccbr.log[ccbr3$minBC > 500,],topoScatter(mean.lphi,sd,log="y",pch=20,resolution=40, 
		xlab=expression(E(log[10](phi))),ylab=expression(sigma),maxFreq=35,thresh=3
	))
	par(op)
},paste0(outdir,"tileSEQ_",geneName,"_errorRegularizationInput"),10,5)
# invisible(dev.off())


logger$info("Performing Error Regularization")

splinemat <- data.frame(
	logsd=log10(ccbr.log$sd),
	logminbc=log10(ccbr3$minBC),
	lphi=ccbr.log$mean.lphi
)
z <- lm(splinemat$logsd ~.,splinemat[,-1])
sdVpred <- 10^cbind(empiric=splinemat$logsd,model=predict(z,splinemat[,-1]))

bnl <- function(pseudo.n,n,model.sd,empiric.sd) {
	sqrt((pseudo.n * model.sd^2 + (n - 1) * empiric.sd^2)/(pseudo.n + n - 2))
}

bayes.sd <- bnl(4,2,sdVpred[,"model"],sdVpred[,"empiric"])

logger$info("Plotting Error Regularization Results")

# pdf(paste0(outdir,"tileSEQ_",geneName,"_errorRegularizationOutput.pdf"),10,5)
html$subsection("Bayesian Regularization of Error")
html$figure(function(){
	op <- par(mfrow=c(1,2))
	topoScatter(sdVpred[,1]+0.0001,sdVpred[,2]+0.0001,log="xy",resolution=50,maxFreq=30,#pch=20,
		xlab=expression("Model"~sigma),ylab=expression("Empiric"~sigma)
	)
	# plot(ccbr.log$sd,bayes.sd,pch=16,col=rgb(79,148,205,50,maxColorValue=255),xlim=c(0,.8),ylim=c(0,.8))
	topoScatter(ccbr.log$sd+0.0001,bayes.sd+0.0001,resolution=60,maxFreq=30,log="xy",
		xlab=expression("Empiric"~sigma),ylab=expression("Bayesian Regularized"~sigma)
	)
	par(op)
},paste0(outdir,"tileSEQ_",geneName,"_errorRegularizationOutput"),10,5)
# invisible(dev.off())

ccbr.log$bsd <- bayes.sd
ccbr.log$minBC <- ccbr3$minBC


logger$info("Writing output")
outfile <- paste0(outdir,"compl_tileSEQ_results_",geneName,".csv")
write.table(ccbr.log,outfile,sep=",",row.names=FALSE)
html$subsection("Output")
html$link.data(outfile)

# write.table(ccbr.log,"ccbr_regularized_SUMO1.csv",sep=",",row.names=FALSE)
# write.table(ccbr.log,"ccbr_regularized.csv",sep=",",row.names=FALSE)

html$shutdown()

logger$info("Done.")
