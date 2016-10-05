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
testfeat <- testable[,-c(4,5,6)]
testscores <- testable[,"score"]


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

testfeat.reach <- testable.reach[,-c(4,5,6)]
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
testfeat.ala <- testable.ala[,-c(4,5,6)]
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

	testfeat.random <- testable.random[,-c(4,5,6)]
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
		random.rmsds,breaks=seq(0.43,0.62,0.005),freq=FALSE,
		border=FALSE,col="gray",
		xlab="RMSD",main="SNP-accessible mutations\nyield worse prediction performance"
	)
	abline(v=rmsds,col="firebrick3",lty="dashed")
	p <- pnorm(q=rmsds[[2]],mean=mean(random.rmsds),sd=sd(random.rmsds),lower.tail=FALSE)
	text(rmsds[[1]],20,sprintf("Full POPCode:\nRMSD = %.02f",rmsds[[1]]),pos=4,cex=0.7)
	text(rmsds[[2]],20,sprintf("SNP-accessible:\nRMSD = %.02f\nP = %.04f",rmsds[[2]],p),pos=4,cex=0.7)
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



