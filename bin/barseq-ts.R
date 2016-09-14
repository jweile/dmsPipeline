###################################
# Process raw barcode reads from  #
# Complementation Time series     #
###################################

options(stringsAsFactors=FALSE)

source("lib/resultfile.R")
source("lib/libyogitools.R")
source("lib/cliargs.R")
source("lib/topoScatter.R")
source("lib/liblogging.R")

#Get output directory
outdir <- getArg("outdir",default="workspace/test/")

#Init logger
logger <- new.logger(paste0(outdir,"barseq-ts.log"))

#Set resultfile
html <- new.resultfile(paste0(outdir,"results.html"))
html$section("UBE2I Complementation Timeseries BarSEQ Data Analysis")

##############
# LOAD INPUT #
##############
logger$info("Reading input")

#Import table linking sample ids to corresponding timepoints and conditions
sample.table <- read.delim("res/muxtags_compTS_samples.tsv",stringsAsFactors=FALSE)

#Read the raw barcode count data
data.mat <- read.csv("input/raw_counts_compTS_joint.csv")

#transform into relative barcode frequencies
data.mat2 <- data.mat[-1,sample.table$sample]
data.rel <- apply(data.mat2,2,function(x)x/sum(x))

#Read the total number of cells in each condition from OD measurements
total.cells <- read.delim("input/total_cells_compTS_joint.tsv")
cells.per.sample <- do.call(c,lapply(unlist(do.call(c,lapply(1:nrow(total.cells),function(i)total.cells[i,])))[-1],rep,3))

#calculate absolute number of cells in each strain/sample
data.abs <- t(apply(data.rel,1,function(x)x*cells.per.sample))

logger$info("Pre-processing kiloseq table")
#Load clone table and determine final genotypes
clone.table <- read.csv("res/clones.csv",stringsAsFactors=FALSE)
rownames(clone.table) <- clone.table$id
final.call <- clone.table$aa.calls
names(final.call) <- clone.table$id

#Process long deletion data to correct genotype calls.
delpos <- do.call(rbind,lapply(strsplit(clone.table$deletions,"-"),function(x) if (length(x)==1&&is.na(x)) c(NA,NA) else as.numeric(x)))
foo <- delpos[,1] < 82 & delpos[,2] > 545
foo[is.na(foo)] <- FALSE
null.clones <- clone.table$id[foo]
foo <- delpos[,1] < delpos[,2]
foo[is.na(foo)] <- FALSE
longdel.clones <- clone.table$id[foo]
wt.clones <- clone.table$id[clone.table$dest.plate=="WT"]
final.call[longdel.clones] <- "longdel"
final.call[null.clones] <- "null"

#################
#Calculate CNFAs#
#################
logger$info("Calculating cumulative normalized fitness advantage for each strain")
cnfas <- do.call(cbind,lapply(1:3,function(repl) {

	perm.rep <- data.frame(
		`t0`=data.abs[,which(with(sample.table, timepoint==0 & replicate==repl))],
		`t6`=data.abs[,which(with(sample.table, timepoint==6 & temp==25 & replicate==repl))],
		`t12`=data.abs[,which(with(sample.table, timepoint==12 & temp==25 & replicate==repl))],
		`t24`=data.abs[,which(with(sample.table, timepoint==24 & temp==25 & replicate==repl))],
		`t48`=data.abs[,which(with(sample.table, timepoint==48 & temp==25 & replicate==repl))]
	)

	sel.rep <- data.frame(
		`t0`=data.abs[,which(with(sample.table, timepoint==0 & replicate==repl))],
		`t6`=data.abs[,which(with(sample.table, timepoint==6 & temp==37 & replicate==repl))],
		`t12`=data.abs[,which(with(sample.table, timepoint==12 & temp==37 & replicate==repl))],
		`t24`=data.abs[,which(with(sample.table, timepoint==24 & temp==37 & replicate==repl))],
		`t48`=data.abs[,which(with(sample.table, timepoint==48 & temp==37 & replicate==repl))]
	)

	all.rep <- cbind(perm=perm.rep,sel=sel.rep)

	all.data <- do.call(rbind,lapply(clone.table$id, function(id) {
		cbind(id=id,call=final.call[[id]],all.rep[id,])
	}))

	#Compute hourly growth rates
	t <- unique(sample.table$timepoint)
	hrates <- function(x) sapply(2:5, function(i) (x[i]/x[i-1])^(1/(t[i]-t[i-1])))
	rates <- do.call(rbind,lapply(1:nrow(all.data), function(i) {
		perm <- unlist(all.data[i,paste("perm.t",t,sep="")])
		sel <- unlist(all.data[i,paste("sel.t",t,sep="")])
		c(
			hrates(perm),
			hrates(sel)
		)
	}))
	rownames(rates) <- rownames(all.data)
	#get average growth rates at each timepoint
	av.rates <- apply(rates,2,function(x)mean(fin(x)))
	#compute fitness advantage compared to pool
	phi <- t(apply(rates,1,`/`,av.rates))
	#normalize fitness advantage to permissive temp
	phi.norm <- phi[,5:8]/phi[,1:4]
	#cumulative normalized fitness advantage
	cnfa <- apply(phi.norm,1,prod)
	return(cnfa)	

}))

#Create a big overview data table linking genotype to phenotype
all.data <- data.frame(id=clone.table$id,mut=final.call,cnfa=cnfas,mean.cnfa=apply(cnfas,1,mean),esd=apply(cnfas,1,sd))
all.data$minCount <- apply(data.abs[all.data$id,which(with(sample.table, timepoint==0))],1,min)

##################
#REGULARIZE STDEV#
##################
logger$info("Regularizing standard deviation")

logger$info("--> Visualizing input")
# pdf(paste0(outdir,"compTS_errorRegularizationInput.pdf"),10,5)
html$subsection("Modeling Error Prior")
html$figure(function(){
	op <- par(mfrow=c(1,2))
	# with(all.data,plot(mean.cnfa,esd,log="y"))
	with(all.data,topoScatter(mean.cnfa,esd,log="y",resolution=50,
		xlim=c(0.25,1.7),ylim=c(1e-4,0.7),
		xlab="Cumulative fitness advantage",ylab="Standard deviation"))
	# with(all.data,plot(minCount,esd,log="xy"))
	with(all.data,topoScatter(minCount,esd,log="xy",resolution=50,
		xlab="permissive Barcode Count",ylab="Standard deviation"))
	par(op)
},paste0(outdir,"compTS_errorRegularizationInput"),10,5)
# invisible(dev.off())

logger$info("--> Calculating")
splinemat <- data.frame(
	logsd=log10(all.data$esd),
	logmin=log10(all.data$minCount),
	cnfa=all.data$mean.cnfa
)
z <- lm(splinemat$logsd ~.,splinemat[,-1])
sdVpred <- 10^cbind(empiric=splinemat$logsd,model=predict(z,splinemat[,-1]))
#Baldi&Long
bnl <- function(pseudo.n,n,model.sd,empiric.sd) {
	sqrt((pseudo.n * model.sd^2 + (n - 1) * empiric.sd^2)/(pseudo.n + n - 2))
}
all.data$bsd <- bnl(6,3,sdVpred[,"model"],sdVpred[,"empiric"])

logger$info("--> Visualizing output")
# pdf(paste0(outdir,"compTS_errorRegularizationOutput.pdf"),5,5)
html$subsection("Bayesian Regularization of Error")
html$figure(function(){
	with(all.data,plot(esd,bsd,
		pch=16,col=rgb(79,148,205,50,maxColorValue=255),
		xlim=c(0,.8),ylim=c(0,.8),
		xlab="Empirical Standard Deviation",
		ylab="Regularized Standard Deviation"
	))
},paste0(outdir,"compTS_errorRegularizationOutput"),5,5)
# invisible(dev.off())

#####################
# QUALITY FILTERING #
#####################
logger$info("Filtering out unreliable clones")

#Filter out broken and badly measured clones
good.data <- all.data[all.data$mut != "longdel",]
good.data <- good.data[good.data$minCount > 500,]
good.data <- good.data[good.data$bsd < 0.05,]
good.clones <- with(clone.table,id[freq > 0.6])
good.data <- good.data[good.data$id %in% good.clones,]

#################
# SCORE SCALING #
#################
logger$info("Adjusting fitness scores to null-wt scale")

good.wt.clones <- intersect(wt.clones,good.data$id)
good.null.clones <- intersect(null.clones,good.data$id)

#mean of null clones
mean.null <- mean(good.data[good.null.clones,"mean.cnfa"])
# score on the scale from null to wt
nw.scale <- (mean(good.data[good.wt.clones,"mean.cnfa"]) - mean(good.data[good.null.clones,"mean.cnfa"]))
score <- ((good.data$mean.cnfa - mean.null) / nw.scale)
good.data$score <- score
good.data$score.bsd <- good.data$bsd / nw.scale

#########################
# OUTPUT PER-CLONE DATA #
#########################
logger$info("Writing per-clone output")

outfile <- paste0(outdir,"compl_timeseries_results_byClone.csv")
write.table(good.data,outfile,sep=",",row.names=FALSE)
html$subsection("Per-clone scores")
html$link.data(outfile)



##############################
# CHECK REPLICATE AGREEMENTS #
##############################
logger$info("Checking replicate agreement")

#Replication plots
techrep.pairs <- do.call(rbind,lapply(1:nrow(good.data),function(i) {
	reps <- ((good.data[i,c("cnfa.1","cnfa.2","cnfa.3")]-mean.null)/nw.scale)
	t(combn(reps,2))
}))
techrep.pairs <- apply(techrep.pairs,2,unlist)

#Check agreement of biological replicates
biorep.pairs <- do.call(rbind,tapply(good.data$score,good.data$mut,function(s) {
	if (length(s) > 1) t(combn(s,2)) else NULL
},simplify=FALSE))

# pdf(paste0(outdir,"compTS_replication.pdf"),10,5)
html$subsection("Replicate correlation")
html$figure(function(){
	op <- par(mfrow=c(1,2))
	topoScatter(
		techrep.pairs[,1], techrep.pairs[,2],resolution=60, 
		xlab="Replicate 1 fitness", ylab="Replicate 2 fitness",
		main="Technical replication"
	)
	text(1,2,sprintf("R = %.02f",cor(techrep.pairs)[1,2]),srt=45)

	topoScatter(biorep.pairs[,1],biorep.pairs[,2],resolution=60,
		thresh=3,maxFreq=50,
		xlab="Replicate 1 fitness", ylab="Replicate 2 fitness",
		main="Biological replication"
	)
	text(.7,2,sprintf("R = %.02f",cor(biorep.pairs)[1,2]),srt=45)
	par(op)
},paste0(outdir,"compTS_replication"),10,5)
# invisible(dev.off())

#################################
# CALCULATE PER-MUTATION SCORES #
#################################
logger$info("Calculating mutation-wise scores")

join.datapoints <- function(ms,sds) {
	ws <- (1/sds)/sum(1/sds)
	mj <- sum(ws*ms)
	vj <- sum(ws*(sds^2+ms^2)) -mj^2
	c(mj=mj,sj=sqrt(vj))
}

#Collapse replicated mutations
mtable <- as.df(tapply(1:nrow(good.data),good.data$mut,function(is) {
	if (length(is) == 1) {
		return(with(good.data,
			list(mut=mut[is],score=score[is],sd=score.bsd[is],se=score.bsd[is]/sqrt(3))
		))
	}
	scores <- good.data$score[is]
	sds <- good.data$score.bsd[is]

	nas <- is.na(scores) | is.na(sds)
	scores <- scores[!nas]
	sds <- sds[!nas]

	j <- join.datapoints(scores,sds)
	se <- j[["sj"]]/sqrt(3*length(scores))

	return(with(good.data,
		list(mut=mut[is[[1]]],score=j[["mj"]],sd=j[["sj"]],se=se)
	))
}))

logger$info("Writing mutation-wise output")
outfile <- paste0(outdir,"compl_timeseries_results_byMut.csv")
write.table(mtable,outfile,sep=",",row.names=FALSE)
html$subsection("Per-mutation scores")
html$link.data(outfile)

html$shutdown()

logger$info("Done.")
