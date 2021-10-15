###########################################################
# Performs subsampling analysis of imputation performance #
# based on SNP-accessible and alanine-scanning mutations  #
###########################################################

source("lib/resultfile.R")
source("lib/libyogitools.R")
source("lib/liblogging.R")
source("lib/cliargs.R")
library("randomForest")
library("parallel")

options(stringsAsFactors=FALSE)

#get output directory
outdir <- getArg("outdir",default="workspace/test/")

#Initialize logger
logger <- new.logger(paste0(outdir,"subsampling.log"))

#Set resultfile
html <- new.resultfile(paste0(outdir,"results.html"))
html$section("Simulate SNP-accessible mutations and Alanine scanning")


scinotExpr <- function(x, digits=2) {
    sign <- ""
    if (x < 0) {
        sign <- "-"
        x <- -x
    }
    exponent <- floor(log10(x))
    if (exponent) {
        xx <- round(x / 10^exponent, digits=digits)
        e <- paste(" %*% 10^", as.integer(exponent), sep="")
    } else {
        xx <- round(x, digits=digits)
        e <- ""
    }
    parse(text=paste("P == ",sign, xx, e, sep=""))
}

##############
# LOAD INPUT #
##############
logger$info("Reading input")

#import feature table and format for random forest
featable <- read.csv(paste0(outdir,"UBE2I_featureMatrix.csv"))
cialevels <- c("no","bb","res")
aas <- c("A","V","L","I","M","F","Y","W","R","H","K","D","E","S","T","N","Q","G","C","P")
featable$mut.aa <- factor(featable$mut.aa,levels=aas)
featable$wt.aa <- factor(featable$wt.aa,levels=aas)
featable$wt.aroVali <- factor(featable$wt.aroVali)
featable$mut.aroVali <- factor(featable$mut.aroVali)
featable$minus2.aroVali <- factor(featable$minus2.aroVali)
featable$minus1.aroVali <- factor(featable$minus1.aroVali)
featable$plus1.aroVali <- factor(featable$plus1.aroVali)
featable$plus2.aroVali <- factor(featable$plus2.aroVali)
featable$sec.struc <- factor(featable$sec.struc)
featable$cia.psiKxE <- factor(featable$cia.psiKxE,levels=cialevels)
featable$cia.rangap <- factor(featable$cia.rangap,levels=cialevels)
featable$cia.sumo <- factor(featable$cia.sumo,levels=cialevels)
featable$cia.ranbp2 <- factor(featable$cia.ranbp2,levels=cialevels)
featable$cia.homodimer <- factor(featable$cia.homodimer,levels=cialevels)
featable$cia.sumo.nc <- factor(featable$cia.sumo.nc,levels=cialevels)
featable$cia.e1 <- factor(featable$cia.e1,levels=cialevels)

testable <- featable[!is.na(featable$score),]
testfeat <- testable[,-c(4,5,6,7)]
testscores <- testable[,"score"]

##############################
# PREPARE REACHABILITY INFORMATION
################################

orf <- scan("res/UBE2I_ORF.fa",what="character")[[2]]

logger$info("Calculating SNP-accessible amino acids")

#translation table
new.trt <- function() {
	library("hash")
	ct <- read.delim("res/codontable.txt",header=FALSE)
	aas <- ct[,2]
	codons <- strsplit(ct[,3],"\\|")
	trtable <- hash()
	for (i in 1:nrow(ct)) {
		for (codon in codons[[i]]){
			trtable[[codon]] <- aas[[i]]
		}
	}
	trtable
}
trt <- new.trt()

codons <- sapply(seq(1,nchar(orf),3),function(i)substr(orf,i,i+2))

aaseq <- sapply(codons,function(codon)trt[[codon]])

#calculate all reachable amino acids
reachable <- lapply(codons, function(codon) {
	unique(do.call(c,lapply(1:3,function(i){
		sapply(c("A","C","T","G"), function(n) {
			mcod <- codon
			substr(mcod,i,i) <- n
			trt[[mcod]]
		})
	})))
})
#remove stop codons
for (i in 1:length(reachable)) {
	reachable[[i]] <- setdiff(reachable[[i]],c(aaseq[[i]],"*"))
}
#format as mutation descriptions
reachable.muts <- do.call(c,lapply(1:length(aaseq),function(pos){
	paste0(aaseq[[pos]],pos,reachable[[pos]])
}))
# z <- randomForest(testfeat,testscores,importance=TRUE)


##########
# Test RMSD of regressed Polyphen2
##########
pp2.raw <- read.delim("res/UBE2I_polyphen2.txt",comment.char="#")
rownames(pp2.raw) <- with(pp2.raw,paste0(o_aa1,o_pos,o_aa2))

testable.pp2 <- data.frame(row.names=testable$mut,dms=testable$score,pp2=pp2.raw[testable$mut,"pph2_prob"])
testable.pp2 <- na.omit(testable.pp2)
testable.pp2$inv <- (1-testable.pp2$pp2+0.0001)/1.0002

logit <- function(p) log(p/(1-p))
testable.pp2$logitInv <- logit(testable.pp2$inv)
lmz <- lm(testable.pp2$dms ~ testable.pp2$logitInv)
testable.pp2$regr <- predict(lmz,testable.pp2)

rmsds.pp2 <- numeric()
rmsds.pp2[["regr"]] <- with(testable.pp2,sqrt(mean((dms-regr)^2)))


####################################
#Test decreasing sizes of training sets

testable <- featable[!is.na(featable$score),]
testfeat <- testable[,-c(4,5,6,7)]
testscores <- testable[,"score"]

reachable.testable.rows <- which(testable$mut %in% reachable.muts)
z <- randomForest(testfeat[reachable.testable.rows,],testscores[reachable.testable.rows])
pred <- predict(z,testfeat[-reachable.testable.rows,])
real <- testable[-reachable.testable.rows,"score"]
reachable.rmsd <- sqrt(mean((pred-real)^2))
reachable.completeness <- length(reachable.testable.rows)/nrow(featable)


#matrix completenesses to simulate
completenesses <- seq(0.05, 1, 0.05)
#restrict completenesses to those that leave at least 1/10 for cross-validation
completenesses <- completenesses[which(completenesses < (nrow(testable)-floor(nrow(testable)/10))/nrow(featable))]
#training set sample sizes corresponding to these completenesses
samplesizes <- floor(completenesses * nrow(featable))

rmsds.ss <- do.call(rbind,lapply(samplesizes, function(samplesize) {
	replicate(3, {

		#Random withholding, not by-the-book cross-validation!!
		toRemove <- nrow(testable) - samplesize

		#number of crossvalidation rounds necessary to approx. cover all testable rows
		nrounds <- round(nrow(testable)/(nrow(testable)-samplesize))

		cat("Samplesize:",samplesize,"\tRounds:",nrounds,"\n")

		pb <- txtProgressBar(max=nrounds,style=3); pr <- 0
		sqds <- do.call(c,mclapply(1:nrounds, function(i) {
			withheldLines <- sample(1:nrow(testable),toRemove,replace=FALSE)
			z <- randomForest(testfeat[-withheldLines,],testscores[-withheldLines])
			setTxtProgressBar(pb,pr <<- pr+1)
			pred <- predict(z,testfeat[withheldLines,])
			real <- testable[withheldLines,"score"]
			(pred-real)^2
		},mc.cores=6))
		close(pb)

		sqrt(mean(sqds))
	})
}))
rownames(rmsds.ss) <- completenesses
colnames(rmsds.ss) <- c("r1","r2","r3")

rmsds.ss <- cbind(
	rmsds.ss,
	mean=apply(rmsds.ss,1,mean),
	se=apply(rmsds.ss,1,sd)/sqrt(3)
)

#For sample size 0, use the score mean for every prediction
random.guessing <- sqrt(mean((testable$score-mean(testable$score))^2))

html$subsection("Sample size vs RMSD")
html$figure(function(){
	# pdf("samplesizeVrmsd.pdf",5,5)
	complPercent <- completenesses*100
	plot(
		complPercent,rmsds.ss[,"mean"],xlim=c(0,100),ylim=c(0.3,0.5),pch=20,
		xlab="Map completeness (%)",ylab="RMSD",col="steelblue3"
	)
	points(0,random.guessing,col="steelblue3",pch=20)
	text(0,random.guessing,"naive guessing",pos=4,cex=0.7)
	points(reachable.completeness*100,reachable.rmsd,col="orange",pch=20)
	text(reachable.completeness*100,reachable.rmsd,"SNP-accessible only",pos=4,cex=0.7)
	arrows(
		complPercent,rmsds.ss[,"mean"]-rmsds.ss[,"se"],
		complPercent,rmsds.ss[,"mean"]+rmsds.ss[,"se"],
		length=0.01,angle=90,code=3,
		col="steelblue3"
	)
	abline(h=rmsds.pp2[[1]],col="orange",lty="dotted")
	text(0,rmsds.pp2[[1]],"PolyPhen-2",pos=4,cex=0.7)

	lfit <- lm(rmsds.ss[,"mean"] ~ complPercent)
	coef <- coefficients(lfit)

	lines(
		seq(0,100,10),coef[[1]]+coef[[2]]*seq(0,100,10),
		col="gray",lty="dashed"
	)
	# dev.off()
},"samplesizeVrmsd",5,5)


####################################
#CROSS-VALIDATION OF COMPLETE MATRIX
logger$info("Replicating 10x cross-validation of random forest imputation")

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

rmsds <- numeric()
rmsds[["POPCode"]] <- sqrt(mean(sqds$sqd))


logger$info("Testing performance of SNP-accessible mutations")
logger$info("-> Recalibrating postional averages")

#TEST MATRIX BASED ON SNP-REACHABLE MUTATIONS
reach.is <- which(testable$mut %in% reachable.muts)
testable.reach <- testable[,-(9:12)]
#re-calculate the positional average based on the smaller subset
testable.reach$posAverage <- sapply(1:nrow(testable.reach),function(k) {
	pos <- testable.reach$pos[[k]]
	mut <- testable.reach$mut[[k]]
	j <- which(testable.reach$mut==mut)
	is <- setdiff(which(testable.reach$pos == pos),j)
	is <- intersect(reach.is,is)
	if (length(is)==0) return(NA)
	else if (length(is)==1) {
		testable.reach$score[is]
	} else {
		scores <- testable.reach$score[is]
		sds <- testable.reach$sd[is]
		weights <- (1/sds)/sum(1/sds)
		sum(scores * weights)
	}
})
testable.reach$posAverageAvail <- !is.na(testable.reach$posAverage)
testable.reach$posAverage[is.na(testable.reach$posAverage)] <- mean(testable.reach$posAverage,na.rm=TRUE)

logger$info("-> Imputing")

testfeat.reach <- testable.reach[,-c(4,5,6,7)]
testscores.reach <- testable.reach[,"score"]

z.reach <- randomForest(testfeat.reach[reach.is,-(6:9)],testscores.reach[reach.is])
pred <- predict(z.reach,testfeat.reach[-reach.is,-(6:9)])
real <- testscores.reach[-reach.is]
# plot(real,pred)
rmsds[["SNP-reachable"]] <- sqrt(mean((pred-real)^2))


###########################
# SIMULATE ALANINE SCANNING
###########################

logger$info("Testing performance of Alanine scan")

ala.is <- which(testable$mut.aa == "A")
testable.ala <- testable[,-(9:12)]
logger$info("-> Recalibrating postional averages")
testable.ala$posAverage <- sapply(1:nrow(testable.ala),function(k) {
	pos <- testable.ala$pos[[k]]
	mut <- testable.ala$mut[[k]]
	j <- which(testable.ala$mut==mut)
	is <- setdiff(which(testable.ala$pos == pos),j)
	is <- intersect(ala.is,is)
	if (length(is)==0) return(NA)
	else if (length(is)==1) {
		testable.ala$score[is]
	} else {
		scores <- testable.ala$score[is]
		sds <- testable.ala$sd[is]
		weights <- (1/sds)/sum(1/sds)
		sum(scores * weights)
	}
})
testable.ala$posAverageAvail <- !is.na(testable.ala$posAverage)
testable.ala$posAverage[is.na(testable.ala$posAverage)] <- mean(testable.ala$posAverage,na.rm=TRUE)

logger$info("-> Imputing")
testfeat.ala <- testable.ala[,-c(4,5,6,7)]
testscores.ala <- testable.ala[,"score"]
z.ala <- randomForest(testfeat.ala[ala.is,-(6:9)],testscores.ala[ala.is])
pred <- predict(z.ala,testfeat.ala[-ala.is,-(6:9)])
real <- testscores.ala[-ala.is]
# plot(real,pred)
rmsds[["Ala-Scan"]] <- sqrt(mean((pred-real)^2))



###########################
# SUBSAMPLING
###########################

logger$info("Testing performance of 500 random subsamples")

#Run 500 mini experiments, drawing random sub-samples of popcode mutations
random.rmsds <- mclapply(1:500,function(x) {
	cat(x,", ")
	testable.random <- testable[,-(9:12)]
	random.is <- sample(1:nrow(testable.random),length(reachable.muts))
	#re-calculate the positional averages based on the random subset
	testable.random$posAverage <- sapply(1:nrow(testable.random),function(k) {
		pos <- testable.random$pos[[k]]
		mut <- testable.random$mut[[k]]
		j <- which(testable.random$mut==mut)
		is <- setdiff(which(testable.random$pos == pos),j)
		is <- intersect(random.is,is)
		if (length(is)==0) return(NA)
		else if (length(is)==1) {
			testable.random$score[is]
		} else {
			scores <- testable.random$score[is]
			sds <- testable.random$sd[is]
			weights <- (1/sds)/sum(1/sds)
			sum(scores * weights)
		}
	})
	testable.random$posAverageAvail <- !is.na(testable.random$posAverage)
	testable.random$posAverage[is.na(testable.random$posAverage)] <- mean(testable.random$posAverage,na.rm=TRUE)

	testfeat.random <- testable.random[,-c(4,5,6,7)]
	testscores.random <- testable.random[,"score"]

	z.random <- randomForest(testfeat.random[random.is,-(6:9)],testscores.random[random.is])
	pred <- predict(z.random,testfeat.random[-random.is,-(6:9)])
	real <- testscores.random[-random.is]
	sqrt(mean((pred-real)^2))
},mc.cores=6)
random.rmsds <- do.call(c,random.rmsds)
cat("\n")




# pdf("subsampling/snpVsPop.pdf",7,4)
html$subsection("SNP-accessible and Alanine Scanning vs POPCode")
html$figure(function(){
	hist(
		random.rmsds,breaks=seq(0.3,0.6,0.005),freq=FALSE,
		border=FALSE,col="gray",
		xlab="RMSD",main="SNP-accessible mutations\nyield worse prediction performance"
	)
	abline(v=rmsds,col="firebrick3",lty="dashed")
	p <- pnorm(q=rmsds[[2]],mean=mean(random.rmsds),sd=sd(random.rmsds),lower.tail=FALSE)
	text(rmsds[[1]],20,sprintf("Full POPCode:\nRMSD = %.02f",rmsds[[1]]),pos=4,cex=0.7)
	text(rmsds[[2]],20,sprintf("SNP-accessible:\nRMSD = %.02f",rmsds[[2]]),pos=4,cex=0.7)
	text(rmsds[[2]],16,scinotExpr(p),pos=4,cex=0.7)

	text(rmsds[[3]],20,sprintf("Alanine-Scan:\nRMSD = %.02f",rmsds[[3]]),pos=4,cex=0.7)
	text(mean(random.rmsds),5,"Random POPCode mut.\nat SNP-accessible sample size",cex=0.7)
},paste0(outdir,"snpVpop"),7,4)
# dev.off()


# #How does alanine scan compare to median position fitness
# median.vals <- with(featable,tapply(score,pos,median,na.rm=TRUE))
# ala.vals <- with(featable, score[mut.aa=="A"])
# layout(cbind(1,2),widths=c(8,1))
# op <- par(mar=c(5,4,4,0)+.1)
# xs <- barplot(ala.vals-median.vals,beside=TRUE,
# 	xlab="AA position",ylab="Ala. Fitness - Median Fitness",
# 	col="steelblue3",border=NA,axes=FALSE,names.arg=NA,space=0
# )
# ticks <- c(1,seq(10,159,10))
# axis(1,at=xs[ticks,],labels=ticks)
# axis(2)
# grid(NA,NULL)
# par(op)
# hdat <- hist(ala.vals-median.vals,breaks=20,plot=FALSE)
# op <- par(mar=c(5,0,4,1)+.1)
# barplot(hdat$density,horiz=TRUE,space=0,
# 	col="steelblue2",border=NA,
# 	xlab="Density",main=""
# )
# par(op)


median.vals <- with(featable,tapply(score,pos,median,na.rm=TRUE))
ala.vals <- with(featable, score[mut.aa=="A"])
ala.diff <- mean(ala.vals-median.vals,na.rm=TRUE)
random.diffs <- replicate(1000,{
	random.vals <- with(featable,tapply(score,pos,function(x)sample(na.omit(x),1)))
	mean(random.vals-median.vals,na.rm=TRUE)
})

# pdf()
html$subsection("SNP-accessible and Alanine Scanning vs POPCode")
html$figure(function(){
	hist(random.diffs,breaks=seq(-0.1,0.3,0.01),col="gray",border=FALSE,freq=FALSE,
		xlab="single fitness - median fitness",main=""
	)
	abline(v=ala.diff,col="firebrick3",lty="dashed")
	p <- pnorm(q=ala.diff,mean=mean(random.diffs),sd=sd(random.diffs),lower.tail=FALSE)
	text(ala.diff,5,sprintf("Alanine Scan:\nFitness difference = %.02f\nP = %.04f",ala.diff,p),pos=4,cex=0.7)
	text(mean(random.diffs),2,"Random Single AA Changes",cex=0.7)
},paste0(outdir,"alaScan"),7,4)




########################################################################################
# Regularize/Impute TileSEQ data and evaluate against HiQual BarSEQ SNP-reachable muts
#  -> Train only on TileSEQ SNPs
#  -> Train on all TileSEQ data
###########################################


#Import TileSEQ and BarSEQ data separately
ccbr <- read.csv(paste0(outdir,"compl_tileSEQ_results_UBE2I_transformed.csv"))
rownames(ccbr) <- ccbr$mut
ltri <- read.csv(paste0(outdir,"compl_timeseries_results_byMut.csv"))

#Define as gold standard those BarSEQ variants that are very well-measured
goldStandard <- ltri[with(ltri,se < 0.07 & regexpr(",",mut) < 1 & !(mut %in% c("null","WT"))),]
goldStandard$pos <- with(goldStandard,as.numeric(substr(mut,2,nchar(mut)-1)))
goldStandard$mutAA <- with(goldStandard,substr(mut,nchar(mut),nchar(mut)))

#gold standard for SNP-accessible mutations
goldenSNPs <- goldStandard[sapply(1:nrow(goldStandard),function(i) {
	pos <- goldStandard[i,"pos"]
	mutAA <- goldStandard[i,"mutAA"]
	mutAA %in% reachable[[pos]]
}),]
rownames(goldenSNPs) <- goldenSNPs$mut


testable.snp <- featable[,-c(10:13)]
rownames(testable.snp) <- featable$mut

#Replace scores with RegSEQ scores
testable.snp$score <- sapply(ccbr[featable$mut,"score.m"], function(s) if(is.na(s)) NA else if(is.infinite(s)) NA else s)
testable.snp$sd <- ccbr[featable$mut,"score.sd"]*5

#Update positional averages based on TileSEQ SNPs
testable.snp$posAverage <- sapply(1:nrow(testable.snp),function(k) {
	pos <- testable.snp$pos[[k]]
	mut <- testable.snp$mut[[k]]
	is <- setdiff(which(testable.snp$pos == pos & 
		testable.snp$mut %in% reachable.muts & 
		!is.na(testable.snp$score)
	),k)
	# is <- is[!is.na(testable.snp$score[is])]
	if (length(is)==0) return(NA)
	else if (length(is)==1) {
		testable.snp$score[is]
	} else {
		scores <- testable.snp$score[is]
		sds <- testable.snp$sd[is]
		weights <- (1/sds)/sum(1/sds)
		sum(scores * weights)
	}
})
testable.snp$posAverageAvail <- !is.na(testable.snp$posAverage)
testable.snp$posAverage[is.na(testable.snp$posAverage)] <- mean(testable.snp$posAverage,na.rm=TRUE)

#cross-validation of goldenSNPs
snps <- goldenSNPs$mut
xvplan <- list()
while(length(snps) > 0) {
	xvplan[[length(xvplan)+1]] <- if (length(snps) >= 25) sample(snps,25) else snps
	snps <- setdiff(snps,xvplan[[length(xvplan)]])
}

#predict snps only based on snp training data
xv <- mclapply(xvplan,function(ms) {
	training <- testable.snp[!is.na(testable.snp$score) & testable.snp$mut %in% reachable.muts,]
	training <- training[-which(training$mut %in% ms),]
	trainfeat <- training[,-c(4,5,6,7)]
	trainscores <- training[,"score"]
	z <- randomForest(trainfeat,trainscores)
	testfeat <- testable.snp[ms,-c(4,5,6,7)]
	pred <- predict(z,testfeat)
	gold <- goldenSNPs[ms,"score"]
	data.frame(mut=ms,gold=gold,ccbr=ccbr[ms,"score.m"],ccbr.sd=ccbr[ms,"score.sd"],pred=pred)
},mc.cores=6)
xv <- do.call(rbind,xv)

#Weighted averages
wa <- function(ms,sds) {
	ws <- (1/sds)/sum(1/sds)
	mj <- sum(ws*ms)
	vj <- sum(ws*(sds^2+ms^2)) -mj^2
	c(mj=mj,sj=sqrt(vj))
}


rf.rmsd <- with(xv, sqrt(mean((gold-pred)^2)))
# ccbr.rmsd <- with(na.omit(xv),sqrt(mean((gold-ccbr)^2)))
regularize <- function(sd2) {
	as.data.frame(do.call(rbind,mapply(function(m1,m2,sd1,sd2){
		if (is.na(m1)) {
			c(mj=m2,sj=sd2)
		} else if (is.na(m2)) {
			c(mj=m1,sj=sd1)
		} else {
			wa(c(m1,m2),c(sd1,sd2))
		}
	},m1=xv$ccbr,m2=xv$pred,sd1=xv$ccbr.sd,sd2=sd2,SIMPLIFY=FALSE)))
}
# regul <-regularize(rf.rmsd)
regul <- regularize(0.06)
regul.rmsd <- sqrt(mean((xv$gold-regul$mj)^2))


#Now predict snps based on all training data

testable.all <- featable[,-c(10:13)]
rownames(testable.all) <- featable$mut

#Replace scores with RegSEQ scores
testable.all$score <- sapply(ccbr[featable$mut,"score.m"], function(s) if(is.na(s)) NA else if(is.infinite(s)) NA else s)
testable.all$sd <- ccbr[featable$mut,"score.sd"]*5

#Update positional averages
testable.all$posAverage <- sapply(1:nrow(testable.all),function(k) {
	pos <- testable.all$pos[[k]]
	mut <- testable.all$mut[[k]]
	is <- setdiff(which(testable.all$pos == pos),k)
	is <- is[!is.na(testable.all$score[is])]
	if (length(is)==0) return(NA)
	else if (length(is)==1) {
		testable.all$score[is]
	} else {
		scores <- testable.all$score[is]
		sds <- testable.all$sd[is]
		weights <- (1/sds)/sum(1/sds)
		sum(scores * weights)
	}
})
testable.all$posAverageAvail <- !is.na(testable.all$posAverage)
testable.all$posAverage[is.na(testable.all$posAverage)] <- mean(testable.all$posAverage,na.rm=TRUE)

xv2 <- mclapply(xvplan,function(ms) {
	training <- testable.all[!is.na(testable.all$score),]
	training <- training[-which(training$mut %in% ms),]
	trainfeat <- training[,-c(4,5,6,7)]
	trainscores <- training[,"score"]
	z <- randomForest(trainfeat,trainscores)
	testfeat <- testable.all[ms,-c(4,5,6,7)]
	pred <- predict(z,testfeat)
	gold <- goldenSNPs[ms,"score"]
	data.frame(mut=ms,gold=gold,ccbr=ccbr[ms,"score.m"],ccbr.sd=ccbr[ms,"score.sd"],pred=pred)
},mc.cores=6)
xv2 <- do.call(rbind,xv2)

rf.rmsd2 <- with(xv2, sqrt(mean((gold-pred)^2)))
# ccbr.rmsd <- with(na.omit(xv2),sqrt(mean((gold-ccbr)^2)))
regularize <- function(sd2) {
	as.data.frame(do.call(rbind,mapply(function(m1,m2,sd1,sd2){
		if (is.na(m1)) {
			c(mj=m2,sj=sd2)
		} else if (is.na(m2)) {
			c(mj=m1,sj=sd1)
		} else {
			wa(c(m1,m2),c(sd1,sd2))
		}
	},m1=xv2$ccbr,m2=xv2$pred,sd1=xv2$ccbr.sd,sd2=sd2,SIMPLIFY=FALSE)))
}
# regul2 <- regularize(rf.rmsd2)
regul2 <- regularize(0.06)
regul2.rmsd <- sqrt(mean((xv2$gold-regul2$mj)^2))

html$subsection("Regularizion based on SNPs vs POPs")
html$figure(function(){
	op <- par(mar=c(7,4,4,1)+1,las=2)
	barplot(
		c(
			impute.snp=rf.rmsd,impute.pop=rf.rmsd2,
			regularize.snp=regul.rmsd,regularize.pop=regul2.rmsd
		),ylab="RMSD",col=paste0("steelblue",c(3,3,4,4)),border="gray"
	)
	grid(NA,NULL)
	par(op)
},paste0(outdir,"snpVpop_regularize"),4,5)
# barplot(c(screen=ccbr.rmsd,predicted=rf.rmsd,regularized=regul.rmsd))
