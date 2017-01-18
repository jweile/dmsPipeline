
source("lib/libyogitools.R")
source("lib/liblogging.R")
source("lib/topoScatter.R")
source("lib/cliargs.R")
source("lib/genophenogram.R")
source("lib/resultfile.R")


library("hash")

options(stringsAsFactors=FALSE)

#get output directory
outdir <- getArg("outdir",default="workspace/test/")
#get input file
infile <- getArg("infile",default="workspace/test/compl_tileSEQ_results_SUMO1_transformed.csv")
# infile <- getArg("infile",default="workspace/test/compl_joint_results_UBE2I.csv")
#get gene name
geneName <- getArg("geneName",default="SUMO1")
# geneName <- getArg("geneName",default="UBE2I")

#Initialize logger
logger <- new.logger(paste0(outdir,"impute_",geneName,".log"))

flipGOF <- as.logical(getArg("flipGOF",default=FALSE))
if (is.na(flipGOF)) {
	err <- "Illegal argument for flipGOF: Must be TRUE or FALSE"
	logger$fatal(err)
	stop(err)
}
#tag for filenames
flptag <- if (flipGOF) "_flipped" else ""

#Set resultfile
html <- new.resultfile(paste0(outdir,"results.html"))
html$section(paste(geneName,"Imputation and regularization",flptag))


#A file that contains a set of mutations that should be excluded from training
# so their imputation results can be evaluated against spotting assays
ctrlSetFile <- getArg("ctrlSet",default=NA)
if (!is.na(ctrlSetFile)) {
	logger$warn("Control set exclusion enabled!")
	ctrlSet <- scan(ctrlSetFile,what="character")
} else {
	ctrlSet <- character()
}

logger$info("Reading input")

screen.data <- read.csv(infile)
rownames(screen.data) <- screen.data$mut

if (geneName == "SUMO1") {

	screen.data <- screen.data[!is.na(screen.data$score.m),]
	screen.data <- screen.data[is.finite(screen.data$score.m),]

	if (flipGOF) {
		flipped <- sapply(screen.data$score.m,function(s) {
			if (s > 1) 1/s else s
		})
		screen.data$score.m <- flipped
	}

} else if (flipGOF) {

	flipped <- sapply(screen.data$score,function(s) {
		if (s > 1) 1/s else s
	})
	screen.data$score <- flipped

}



prot <- scan(paste0("res/",geneName,"_aa.fa"),what="character",sep="\n")[[2]]
wt.aa <- sapply(1:nchar(prot),function(i)substr(prot,i,i))

conservation <- as.integer(scan(paste0("res/",geneName,"_cons.txt"),what="integer",sep=","))
aas <- c("A","V","L","I","M","F","Y","W","R","H","K","D","E","S","T","N","Q","G","C","P")


provean <- read.csv("res/provean.csv")
provean <- provean[provean$protein==geneName,]
rownames(provean) <- provean$mut


###########################
# CONSTRUCT FEATURE TABLE #
###########################

logger$info("Building feature table.")

featable <- expand.grid(pos=1:length(wt.aa),mut.aa=aas)
featable$wt.aa <- wt.aa[featable$pos]
featable$mut <- with(featable,paste(wt.aa,pos,mut.aa,sep=""))

if (geneName == "UBE2I") {
	featable$score <- sapply(featable$mut,function(m) screen.data[m,"score"])
	featable$sd <- sapply(featable$mut,function(m) screen.data[m,"sd"])
} else {
	featable$score <- sapply(featable$mut,function(m) screen.data[m,"score.m"])
	featable$sd <- sapply(featable$mut,function(m) screen.data[m,"score.sd"])
}
featable$score[featable$wt.aa == featable$mut.aa] <- 1
featable$sd[featable$wt.aa == featable$mut.aa] <- 0.0001

screen.singles <- screen.data[!(screen.data$mut %in% c("WT","null")) & regexpr(",",screen.data$mut) < 1,]
screen.singles <- screen.singles[regexpr("_",screen.singles$mut) < 1,]
screen.singles$pos <- as.integer(substr(screen.singles$mut,2,nchar(screen.singles$mut)-1))

if (geneName == "UBE2I") {

	logger$info(" -> Calculating positional averages")

	featable$posAverage <- sapply(1:nrow(featable),function(k) {
		pos <- featable$pos[[k]]
		mut <- featable$mut[[k]]
		j <- which(screen.singles$mut==mut)
		is <- setdiff(which(screen.singles$pos == pos),j)
		if (length(is)==0) return(NA)
		else if (length(is)==1) {
			screen.singles$score[is]
		} else {
			scores <- screen.singles$score[is]
			sds <- screen.singles$sd[is]
			weights <- (1/sds)/sum(1/sds)
			sum(scores * weights)
		}
	})
	featable$posAverageAvail <- !is.na(featable$posAverage)
	featable$posAverage[is.na(featable$posAverage)] <- mean(featable$posAverage,na.rm=TRUE)

	logger$info(" -> Calculating multi-mutant averages")

	multimuts <- strsplit(screen.data$mut,",")
	anymut.idx <- hash()
	for (i in 1:nrow(screen.data)) {
		if (length(multimuts[[i]] > 1)) {
			for (m in multimuts[[i]]) {
				anymut.idx[[m]] <- c(anymut.idx[[m]],i)
			}
		}
	}
	singlemut.idx <- hash()
	for (i in 1:nrow(screen.singles)) {
		m <- screen.singles$mut[[i]]
		singlemut.idx[[m]] <- c(singlemut.idx[[m]],i)
	}

	featable$mmAverage <- sapply(featable$mut, function(m) {
		is <- anymut.idx[[m]]
		if (length(is) > 2) {
			scores <- screen.data$score[is]
			sds <- screen.data$sd[is]
			weights <- (1/sds)/sum(1/sds)
			sum(scores * weights)
		} else NA
	})
	featable$mmAverageAvail <- !is.na(featable$mmAverage)
	featable$mmAverage[is.na(featable$mmAverage)] <- mean(featable$mmAverage,na.rm=TRUE)


	logger$info(" -> Calculating multiplicative model predictions")

	featable$giImpute <- sapply(featable$mut, function(m) {
		is <- anymut.idx[[m]]
		preds <- do.call(c,lapply(is, function(i) {
			dm.score <- screen.data[i,"score"]
			other.ms <- setdiff(multimuts[[i]],m)
			if (length(other.ms) == 1 && has.key(other.ms,singlemut.idx)) {
				sm.score <- screen.singles[singlemut.idx[[other.ms]],"score"]
				sms <- mean(sm.score,na.rm=TRUE)
				if (dm.score <= 0) {
					if (sms <= 0) NA else 0
				} else if (sms <= 0) {
					NA
				} else {
					out <- dm.score / sms
					if (out > 3) 3 else out
				}
			} else NULL
		}))
		mean(preds,na.rm=TRUE)
	})
	featable$giImputeAvail <- !is.na(featable$giImpute)
	featable$giImpute[is.na(featable$giImpute)] <- mean(featable$giImpute,na.rm=TRUE)

} else {

	logger$info(" -> Calculating positional averages")

	screen.singles <- screen.singles[!is.na(screen.singles$score.m),]
	screen.singles <- screen.singles[!is.infinite(screen.singles$score.m),]

	featable$posAverage <- sapply(1:nrow(featable),function(k) {
		pos <- featable$pos[[k]]
		mut <- featable$mut[[k]]
		j <- which(screen.singles$mut==mut)
		is <- setdiff(which(screen.singles$pos == pos),j)
		if (length(is)==0) return(NA)
		else if (length(is)==1) {
			screen.singles$score.m[is]
		} else {
			scores <- screen.singles$score.m[is]
			sds <- screen.singles$score.sd[is]
			weights <- (1/sds)/sum(1/sds)
			sum(scores * weights)
		}
	})
	featable$posAverageAvail <- !is.na(featable$posAverage)
	featable$posAverage[is.na(featable$posAverage)] <- mean(featable$posAverage,na.rm=TRUE)

}

featable$conservation <- conservation[featable$pos]

library("Biostrings")
data(BLOSUM62)
featable$blosum <- mapply(function(from,to){
	BLOSUM62[from,to]
},from=featable$wt.aa,to=featable$mut.aa)



provean.feat <- provean[featable$mut,c("provean","sift")]
if (any(is.na(provean.feat$provean))){
	provean.feat[which(is.na(provean.feat$provean)),"provean"] <- mean(provean.feat$provean,na.rm=TRUE)
}
if (any(is.na(provean.feat$sift))){
	provean.feat[which(is.na(provean.feat$sift)),"sift"] <- mean(provean.feat$sift,na.rm=TRUE)
}
featable <- cbind(featable,provean.feat)


logger$info(" -> Assigning biochemical properties")

aa.props <- read.csv("res/aa_props.csv")
rownames(aa.props) <- aa.props[,"AA"]

featable$wt.mass <- aa.props[featable$wt.aa,"Mass"]
featable$mut.mass <- aa.props[featable$mut.aa,"Mass"]
featable$diff.mass <- with(featable,mut.mass - wt.mass)

featable$wt.volume <- aa.props[featable$wt.aa,"Volume"]
featable$mut.volume <- aa.props[featable$mut.aa,"Volume"]
featable$diff.volume <- with(featable,mut.volume - wt.volume)

featable$wt.isoelec <- aa.props[featable$wt.aa,"Isoelecticity"]
featable$mut.isoelec <- aa.props[featable$mut.aa,"Isoelecticity"]
featable$diff.isoelec <- with(featable,mut.isoelec - wt.isoelec)

featable$wt.pkac <- aa.props[featable$wt.aa,"Carboxyl.dissociation"]
featable$mut.pkac <- aa.props[featable$mut.aa,"Carboxyl.dissociation"]
featable$diff.pkac <- with(featable,mut.pkac - wt.pkac)

featable$wt.pkaa <- aa.props[featable$wt.aa,"Amino.dissociation"]
featable$mut.pkaa <- aa.props[featable$mut.aa,"Amino.dissociation"]
featable$diff.pkaa <- with(featable,mut.pkaa - wt.pkaa)

featable$wt.pkar <- aa.props[featable$wt.aa,"Residue.dissociation"]
featable$mut.pkar <- aa.props[featable$mut.aa,"Residue.dissociation"]
featable$diff.pkar <- with(featable,mut.pkar - wt.pkar)

featable$wt.hydropathy <- aa.props[featable$wt.aa,"Hydropathy"]
featable$mut.hydropathy <- aa.props[featable$mut.aa,"Hydropathy"]
featable$diff.hydropathy <- with(featable,mut.hydropathy - wt.hydropathy)

featable$wt.polar <- aa.props[featable$wt.aa,"Polar"]
featable$mut.polar <- aa.props[featable$mut.aa,"Polar"]
featable$diff.polar <- with(featable,mut.polar != wt.polar)

featable$wt.charge <- aa.props[featable$wt.aa,"Charge"]
featable$mut.charge <- aa.props[featable$mut.aa,"Charge"]
featable$diff.charge <- with(featable,mut.charge - wt.charge)

featable$wt.acidity <- aa.props[featable$wt.aa,"acidic.basic"]
featable$mut.acidity <- aa.props[featable$mut.aa,"acidic.basic"]
featable$diff.acidity <- with(featable,mut.acidity - wt.acidity)

featable$wt.aroVali <- aa.props[featable$wt.aa,"Aromatic.or.Aliphatic"]
featable$mut.aroVali <- aa.props[featable$mut.aa,"Aromatic.or.Aliphatic"]

featable$wt.incrate <- aa.props[featable$wt.aa,"Incorporation.rate"]
featable$mut.incrate <- aa.props[featable$mut.aa,"Incorporation.rate"]
featable$diff.incrate <- with(featable,mut.incrate - wt.incrate)

featable$wt.burialRate <- aa.props[featable$wt.aa,"Burial.rate"]
featable$mut.burialRate <- aa.props[featable$mut.aa,"Burial.rate"]
featable$diff.burialRate <- with(featable,mut.burialRate - wt.burialRate)

featable$wt.avAcc <- aa.props[featable$wt.aa,"Average.Accessibility"]
featable$mut.avAcc <- aa.props[featable$mut.aa,"Average.Accessibility"]
featable$diff.avAcc <- with(featable,mut.avAcc - wt.avAcc)

offsets <- as.df(lapply(featable$pos,function(i) {
	list(
		minus2= if(i <= 2) i+2 else i-2,
		minus1= if(i <= 1) i+1 else i-1,
		plus1 = if(i > length(wt.aa)-1) i-1 else i+1,
		plus2 = if(i > length(wt.aa)-2) i-2 else i+2
	)
}))

foo <- apply(offsets,2,function(i) aa.props[wt.aa[i],"Mass"])
colnames(foo) <- paste0(colnames(foo),".mass")
featable <- cbind(featable,foo)

foo <- apply(offsets,2,function(i) aa.props[wt.aa[i],"Volume"])
colnames(foo) <- paste0(colnames(foo),".volume")
featable <- cbind(featable,foo)

foo <- apply(offsets,2,function(i) aa.props[wt.aa[i],"Isoelecticity"])
colnames(foo) <- paste0(colnames(foo),".isoelec")
featable <- cbind(featable,foo)

foo <- apply(offsets,2,function(i) aa.props[wt.aa[i],"Carboxyl.dissociation"])
colnames(foo) <- paste0(colnames(foo),".pkac")
featable <- cbind(featable,foo)

foo <- apply(offsets,2,function(i) aa.props[wt.aa[i],"Amino.dissociation"])
colnames(foo) <- paste0(colnames(foo),".pkaa")
featable <- cbind(featable,foo)

foo <- apply(offsets,2,function(i) aa.props[wt.aa[i],"Residue.dissociation"])
colnames(foo) <- paste0(colnames(foo),".pkar")
featable <- cbind(featable,foo)

foo <- apply(offsets,2,function(i) aa.props[wt.aa[i],"Hydropathy"])
colnames(foo) <- paste0(colnames(foo),".hydropathy")
featable <- cbind(featable,foo)

foo <- apply(offsets,2,function(i) aa.props[wt.aa[i],"Polar"])
colnames(foo) <- paste0(colnames(foo),".polar")
featable <- cbind(featable,foo)

foo <- apply(offsets,2,function(i) aa.props[wt.aa[i],"Charge"])
colnames(foo) <- paste0(colnames(foo),".charge")
featable <- cbind(featable,foo)

foo <- apply(offsets,2,function(i) aa.props[wt.aa[i],"acidic.basic"])
colnames(foo) <- paste0(colnames(foo),".acidity")
featable <- cbind(featable,foo)

foo <- apply(offsets,2,function(i) aa.props[wt.aa[i],"Aromatic.or.Aliphatic"])
colnames(foo) <- paste0(colnames(foo),".aroVali")
featable <- cbind(featable,foo)

foo <- apply(offsets,2,function(i) aa.props[wt.aa[i],"Incorporation.rate"])
colnames(foo) <- paste0(colnames(foo),".incrate")
featable <- cbind(featable,foo)

foo <- apply(offsets,2,function(i) aa.props[wt.aa[i],"Burial.rate"])
colnames(foo) <- paste0(colnames(foo),".burialRate")
featable <- cbind(featable,foo)

foo <- apply(offsets,2,function(i) aa.props[wt.aa[i],"Average.Accessibility"])
colnames(foo) <- paste0(colnames(foo),".avAcc")
featable <- cbind(featable,foo)


logger$info(" -> Loading structural features")

strucfeats <- read.csv(paste0("res/",geneName,"_structural_features.csv"))

featable <- cbind(featable,strucfeats[featable$pos,2:ncol(strucfeats)])

#Convert strings to factors to accomodate for randomForest input requirements
featable$mut.aa <- factor(featable$mut.aa,levels=aas)
featable$wt.aa <- factor(featable$wt.aa,levels=aas)
featable$wt.aroVali <- factor(featable$wt.aroVali)
featable$mut.aroVali <- factor(featable$mut.aroVali)
featable$minus2.aroVali <- factor(featable$minus2.aroVali)
featable$minus1.aroVali <- factor(featable$minus1.aroVali)
featable$plus1.aroVali <- factor(featable$plus1.aroVali)
featable$plus2.aroVali <- factor(featable$plus2.aroVali)
featable$sec.struc <- factor(featable$sec.struc)

clevels <- c("no","bb","res")

if (geneName == "UBE2I") {

	featable$cia.psiKxE <- factor(featable$cia.psiKxE,levels=clevels)
	featable$cia.rangap <- factor(featable$cia.rangap,levels=clevels)
	featable$cia.sumo <- factor(featable$cia.sumo,levels=clevels)
	featable$cia.ranbp2 <- factor(featable$cia.ranbp2,levels=clevels)
	featable$cia.homodimer <- factor(featable$cia.homodimer,levels=clevels)
	featable$cia.sumo.nc <- factor(featable$cia.sumo.nc,levels=clevels)
	featable$cia.e1 <- factor(featable$cia.e1,levels=clevels)

} else {

	featable$daxx.cia <- factor(featable$daxx.cia,levels=clevels)
	featable$e1.cia <- factor(featable$e1.cia,levels=clevels)
	featable$pias.cia <- factor(featable$pias.cia,levels=clevels)
	featable$pml.cia <- factor(featable$pml.cia,levels=clevels)
	featable$ranbp.cia <- factor(featable$ranbp.cia,levels=clevels)
	featable$senp1.cia <- factor(featable$senp1.cia,levels=clevels)
	featable$senp2.cia <- factor(featable$senp2.cia,levels=clevels)
	featable$tdg.cia <- factor(featable$tdg.cia,levels=clevels)
	featable$ube2i.cov.cia <- factor(featable$ube2i.cov.cia,levels=clevels)
	featable$ubei2i.noncov.cia <- factor(featable$ubei2i.noncov.cia,levels=clevels)
}

featable$numInterfaces <- apply(
	featable[,regexpr("^b\\.",colnames(featable))>0],
	1, function(x) sum(x > 0)
)

featable$numCIA.any <- apply(
	featable[,regexpr("cia",colnames(featable))>0],
	1, function(x) sum(x != "no")
)

featable$numCIA.res <- apply(
	featable[,regexpr("cia",colnames(featable))>0],
	1, function(x) sum(x == "res")
)


logger$info(" -> Saving feature table")
outfile <- paste0(outdir,geneName,flptag,"_featureMatrix.csv")
write.table(featable,outfile,sep=",",row.names=FALSE)
html$subsection("Feature table")
html$link.data(outfile)



############################################
# Baseline naive and regression imputation #
############################################
# logger$info("Calculating Imputation baseline")

testable <- featable[!is.na(featable$score) & !(featable$mut %in% ctrlSet),]


# #Linear regression based on positional average, multimutant and multiplicative numbers

# rmsds <- numeric()
# logger$info(" -> Cross validating regession based on averages")
# rmsds[["averagesOnly"]] <- sqrt(mean(sapply(1:nrow(testable), function(i) {
# 	coef <- coefficients(with(testable[-i,],lm(score ~ posAverage + mmAverage + giImpute)))
# 	pred <- coef[["(Intercept)"]] + sum(sapply(names(coef)[-1],function(k)coef[[k]]*testable[i,k]))
# 	real <- testable[i,"score"]
# 	(pred-real)^2
# })))

# logger$info(" -> Cross validating regession based on averages and BLOSUM")
# rmsds[["blosum"]] <- sqrt(mean(sapply(1:nrow(testable), function(i) {
# 	coef <- coefficients(with(testable[-i,],lm(score ~ posAverage + mmAverage + giImpute + blosum)))
# 	pred <- coef[["(Intercept)"]] + sum(sapply(names(coef)[-1],function(k)coef[[k]]*testable[i,k]))
# 	real <- testable[i,"score"]
# 	(pred-real)^2
# })))

# logger$info(" -> Cross validating regession based on averages and conservation")
# rmsds[["conservation"]] <- sqrt(mean(sapply(1:nrow(testable), function(i) {
# 	coef <- coefficients(with(testable[-i,],lm(score ~ posAverage + mmAverage + giImpute + conservation)))
# 	pred <- coef[["(Intercept)"]] + sum(sapply(names(coef)[-1],function(k)coef[[k]]*testable[i,k]))
# 	real <- testable[i,"score"]
# 	(pred-real)^2
# })))

# logger$info(" -> Cross validating regession based on averages, BLOSUM and conservation")
# rmsds[["blosum+cons"]] <- sqrt(mean(sapply(1:nrow(testable), function(i) {
# 	coef <- coefficients(with(testable[-i,],lm(score ~ posAverage + mmAverage + giImpute + blosum + conservation)))
# 	pred <- coef[["(Intercept)"]] + sum(sapply(names(coef)[-1],function(k)coef[[k]]*testable[i,k]))
# 	real <- testable[i,"score"]
# 	(pred-real)^2
# })))



######################
# PERFORM IMPUTATION #
######################
logger$info("Performing RandomForest Imputation")

library("randomForest")

testfeat <- testable[,-c(4,5,6)]
testscores <- testable[,"score"]
z <- randomForest(testfeat,testscores,importance=TRUE)

# pdf(paste0(outdir,"imputation_",geneName,"_variableImportance.pdf"))
html$subsection("Feature importance")
html$figure(function(){
	varImpPlot(z)
},paste0(outdir,"imputation_",geneName,"_variableImportance",flptag))
# invisible(dev.off())


logger$info(" -> 10x cross validation")
library(parallel)

##Make a 10x cross-validation plan
chunksize <- floor(nrow(testable)/10)
pool <- 1:nrow(testable)
chunks <- replicate(10,{
	chunk <- if (length(pool) > chunksize) sample(pool,chunksize) else pool
	pool <<- setdiff(pool,chunk)
	chunk
},simplify=FALSE)

##Run cross-validation
pb <- txtProgressBar(max=10,style=3); pr <- 0
sqds <- do.call(rbind,mclapply(1:10, function(i) {
	chunk <- chunks[[i]]
	z <- randomForest(testfeat[-chunk,],testscores[-chunk])
	setTxtProgressBar(pb,pr <<- pr+1)
	pred <- predict(z,testfeat[chunk,])
	real <- testable[chunk,"score"]
	pos <- testable[chunk,"pos"]
	mutaa <- testable[chunk,"mut.aa"]
	data.frame(idx=chunk,pos=pos,mut.aa=mutaa,score=real,pred=pred,sqd=(pred-real)^2)
},mc.cores=6))
close(pb)

sqds <- sqds[order(sqds$idx),]
rmsd <- sqrt(mean(sqds$sqd))


# barplot(rmsds)
logger$info(" -> Plotting cross-validation result")

# pdf(paste0(outdir,"imputation_",geneName,"_errorMap.pdf"),6*4,4)
html$subsection("Squared deviation map")
html$figure(function(){
	layout(cbind(1,2),widths=c(9.25,.75))
	#main plot
	op <- par(cex=.6,las=1,mar=c(5,4,0,0)+.1)
	plot(NA,type="n",
		xlim=c(0,length(wt.aa)+1),ylim=c(0,length(aas)+1),axes=FALSE,
		xlab="AA position",ylab="AA residue",main=""
	)
	axis(1,c(1,seq(5,160,5)))
	axis(2,at=1:20,labels=rev(aas))
	x <- sqds$pos
	y <- length(aas) - sapply(sqds$mut.aa,function(a)which(aas==a)) + 1
	colRamp <- colorRampPalette(c("green","red"))(11)
	col.idx <- sapply(round(1+10*sqrt(sqds$sqd)),function(x)min(c(x,11)))
	cols <- colRamp[col.idx]
	# cols <- colRamp[round(10*diffs/max(diffs))+1]
	rect(x-.5,y-.5,x+.5,y+.5,col=cols,border=NA)
	par(op)
	#legend
	op <- par(mar=c(5,0,1,4)+.1)
	plot(NA,type="n",xlim=c(0,1),ylim=c(0,11),axes=FALSE,xlab="",ylab="")
	rect(0,0:10,1,1:11,col=colRamp,border=NA)
	axis(4,at=c(.5,10.5),labels=c("0",">=1"))
	mtext("error",side=4,line=2)
	par(op)
},paste0(outdir,"imputation_",geneName,"_errorMap",flptag),6*4,4)
# invisible(dev.off())

html$subsection("Scores vs RandomForest predictions")
html$figure(function(){
	if (geneName == "UBE2I") {
		# pdf(paste0(outdir,"imputation_",geneName,"_scoreVpred.pdf"),9,4)
		#Plot growth score vs predicted score
		with(sqds,topoScatter(score,pred,xlab="real score",ylab="predicted score",
			main=sprintf("R = %.2f",cor(score,pred)),resolution=60
		))
		abline(v=0:1,h=0:1,col=c("firebrick3","darkolivegreen3"))
		# invisible(dev.off())
	} else {
		# pdf(paste0(outdir,"imputation_",geneName,"_scoreVpred.pdf"),7,4)
		#Plot growth score vs predicted score
		with(sqds,plot(score,pred,xlab="real score",ylab="predicted score",
			main=sprintf("R = %.2f",cor(score,pred)),pch=20
		))
		abline(v=0:1,h=0:1,col=c("firebrick3","darkolivegreen3"))
		# invisible(dev.off())
	}
},paste0(outdir,"imputation_",geneName,"_scoreVpred",flptag),if(geneName=="UBE2I") 9 else 7, 4)


#############################
# Merge data and prediction #
#############################
logger$info("Regularizing and completing data using prediction")

#Build complete joint table
join.datapoints <- function(ms,sds) {
	ws <- (1/sds)/sum(1/sds)
	mj <- sum(ws*ms)
	vj <- sum(ws*(sds^2+ms^2)) -mj^2
	c(mj=mj,sj=sqrt(vj))
}


score.table <- data.frame(mut=featable$mut,screen.score=featable$score,screen.sd=featable$sd)

feat.in <- featable[,-c(4,5,6)]
score <- predict(z,feat.in)
score.table$predicted.score <- score
score.table$joint.score <- score
score.table$joint.sd <- rmsd
for (i in which(!is.na(featable$score))) {
	vs <- c(pred=score[[i]],screen=featable$score[[i]])
	ss <- c(pred=rmsd,screen=featable$sd[[i]])
	joint <- join.datapoints(vs,ss)
	score[[i]] <- joint[["mj"]]
	score.table$joint.score[[i]] <- joint[["mj"]]
	score.table$joint.sd[[i]] <- joint[["sj"]]
}

logger$info("Writing result to file.")

outfile <- paste0(outdir,"imputed_regularized_",geneName,flptag,"_scores.csv")
write.table(score.table,outfile,sep=",",row.names=FALSE)
html$subsection("Output")
html$link.data(outfile)


logger$info("Drawing complete genophenogram")

# pdf(paste0(outdir,"imputed_regularized_",geneName,"_genophenogram.pdf"),19,4)
html$subsection("Regularized and Imputed Genophenogram")
html$figure(function(){
	genophenogram(wt.aa,featable$pos,featable$mut.aa,score.table$joint.score)
},paste0(outdir,"imputed_regularized_",geneName,"_genophenogram",flptag),if(geneName=="UBE2I")19 else 16,4)
# invisible(dev.off())

html$shutdown()

logger$info("Done.")
