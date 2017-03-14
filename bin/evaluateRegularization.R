options(stringsAsFactors=FALSE)

source("lib/resultfile.R")
source("lib/liblogging.R")
source("lib/cliargs.R")
# source("lib/libyogitools.R")

library("randomForest")


outdir <- getArg("outdir",default="workspace/test/")
# outdir <- "workspace/20170131-175115/"

#Initialize logger
logger <- new.logger(paste0(outdir,"evaluateRegularization.log"))

#Set resultfile
html <- new.resultfile(paste0(outdir,"results.html"))
html$section("Evaluate Regularization")

logger$info("Loading data...")

#read barseq data
barseq <- read.csv(paste0(outdir,"compl_timeseries_results_byMut.csv"))
rownames(barseq) <- barseq$mut
#read tileseq data
tileseq <- read.csv(paste0(outdir,"compl_tileSEQ_results_UBE2I_transformed.csv"))
rownames(tileseq) <- tileseq$mut
#remove any invalid values from tileseq
tileseq <- tileseq[is.finite(tileseq$score.m),]
tileseq <- tileseq[regexpr("_",tileseq$mut) < 1,]

#compare minBC with stderr
# with(tileseq,plot(minBC,score.se))
# op <- par(mfrow=c(2,1))
# breaks <- seq(0,0.5,0.005)
# hist(barseq$se,breaks=breaks)
# hist(tileseq$score.se,breaks=breaks)
# abline(v=0.05,col=2)
# par(op)

#find clones that are bad in tileseq but good in barseq
bad.tile.clones <- with(tileseq,mut[which(score.se > 0.05)])
good.bar.clones <- with(barseq,mut[which(se < 0.05)])
test.clones <- intersect(bad.tile.clones,good.bar.clones)

common.clones <- intersect(barseq$mut,tileseq$mut)

# diffs <- abs(barseq[common.clones,"score"]-tileseq[common.clones,"score.m"])
# barse <- log10(barseq[common.clones,"se"])
# tilese <- log10(tileseq[common.clones,"score.se"])
# tilereads <- log10(tileseq[common.clones,"minBC"])
# prod.se <- mapply(`+`,x=barse,y=tilese)
# op <- par(mfrow=c(1,4))
# plot(prod.se,diffs)
# plot(barse,diffs)
# plot(tilese,diffs)
# plot(tilereads,diffs)
# par(op)



featable <- read.csv(paste0(outdir,"UBE2I_featureMatrix.csv"))
rownames(featable) <- featable$mut
clevels <- c("no","bb","res")
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
featable$cia.psiKxE <- factor(featable$cia.psiKxE,levels=clevels)
featable$cia.rangap <- factor(featable$cia.rangap,levels=clevels)
featable$cia.sumo <- factor(featable$cia.sumo,levels=clevels)
featable$cia.ranbp2 <- factor(featable$cia.ranbp2,levels=clevels)
featable$cia.homodimer <- factor(featable$cia.homodimer,levels=clevels)
featable$cia.sumo.nc <- factor(featable$cia.sumo.nc,levels=clevels)
featable$cia.e1 <- factor(featable$cia.e1,levels=clevels)


#Overwrite featable with values from only tileseq
featable$score <- NA
featable$sd <- NA
featable$df <- NA
featable[tileseq$mut,"score"] <- tileseq$score.m
featable[tileseq$mut,"sd"] <- tileseq$score.sd
featable[tileseq$mut,"df"] <- 4 
#restore wildtype values
featable$score[featable$wt.aa == featable$mut.aa] <- 1
featable$sd[featable$wt.aa == featable$mut.aa] <- 0.00001



#Run RandomForest
logger$info("Running RandomForest...")

testable <- featable[!is.na(featable$score),]
testfeat <- testable[,-c(4,5,6,7)]
testscores <- testable[,"score"]
z <- randomForest(testfeat,testscores,importance=FALSE)

logger$info("Regularizing...")

#RMSD achieved by random forest on full dataset
rmsd <- 0.33
# rmsd <- 0.05

#Build complete joint table
join.datapoints <- function(ms,sds,ses) {
	ws <- (1/ses)/sum(1/ses)
	mj <- sum(ws*ms)
	vj <- sum(ws*(sds^2+ms^2)) -mj^2
	c(mj=mj,sj=sqrt(vj))
}

score.table <- data.frame(mut=featable$mut,screen.score=featable$score,
	screen.sd=featable$sd,df=featable$df,row.names=featable$mut)

feat.in <- featable[,-c(4,5,6,7)]
score <- predict(z,feat.in)
score.table$predicted.score <- score
score.table$joint.score <- score
score.table$joint.sd <- rmsd
score.table$joint.se <- rmsd
for (i in which(!is.na(featable$score))) {
	if (featable[i,"mut.aa"]==featable[i,"wt.aa"]) {
		score.table$joint.score[[i]] <- featable$score[[i]]
		score.table$joint.sd[[i]] <- featable$sd[[i]]
		score.table$joint.se[[i]] <- featable$sd[[i]]
	} else {
		vs <- c(pred=score[[i]],screen=featable$score[[i]])
		sds <- c(pred=rmsd,screen=featable$sd[[i]])
		ses <- c(pred=rmsd,screen=featable$sd[[i]]/sqrt(featable$df[[i]]))
		joint <- join.datapoints(vs,sds,ses)
		score[[i]] <- joint[["mj"]]
		score.table$joint.score[[i]] <- joint[["mj"]]
		score.table$joint.sd[[i]] <- joint[["sj"]]
		score.table$joint.se[[i]] <- joint[["sj"]]/sqrt(featable$df[[i]])
	}
}

eval.table <- cbind(score.table[test.clones,],barseq=barseq[test.clones,c("score","se")])

logger$info("Drawing Robin Hood Plot")
html$subsection("Regularization accuracy")
# pdf(paste0(outdir,"regularization_effect.pdf"),3,7)
html$figure(function() {
	n <- length(test.clones)
	op <- par(las=2)
	plot(0,type="n",ylim=c(0,n+1),xlim=c(-0.3,1.2),xlab="score",ylab="variant",axes=FALSE)
	axis(1)
	abline(h=1:n,lty="dotted",col="gray")
	abline(v=0:1,lty="dashed",col=c("firebrick3","darkolivegreen3"))
	text(c(0,1),c(n+1,n+1),c("null","wt"),pos=2,srt=90,col=c("firebrick3","darkolivegreen3"))
	points(eval.table$predicted.score,1:n,col="gray")
	points(eval.table$barseq.score,1:n,col="firebrick3",pch=13)
	with(eval.table,arrows(screen.score,1:n,joint.score,1:n,length=0.04))
	mtext(eval.table$mut,side=2,at=1:n)
	par(op)
	legend("bottomright",
		c("Predicted value","High Quality BarSEQ value","Regularization from/to"),
		col=c("gray","firebrick3","black"),pch=c("\U2B58","\U2B59","\U2192"),cex=0.8
	)
},paste0(outdir,"regularization_effect"),7,4)
# dev.off()


##OLD vertical version of the plot
# n <- length(test.clones)
# op <- par(las=2)
# plot(0,type="n",xlim=c(0,n+1),ylim=c(-0.5,1.2),ylab="score",xlab="variant",axes=FALSE)
# axis(2)
# # grid(NA,NULL)
# abline(v=1:n,lty="dotted",col="gray")
# abline(h=0:1,lty="dashed",col=c("firebrick3","darkolivegreen3"))
# # points(1:n,eval.table$screen.score,col="steelblue3",pch=20)
# points(1:n,eval.table$predicted.score,col="orange")
# # points(1:n,eval.table$joint.score,col="firebrick3",pch=20)
# points(1:n,eval.table$barseq.score,col="chartreuse3",pch=20)
# with(eval.table,arrows(1:n,screen.score,1:n,joint.score,length=0.04))
# mtext(eval.table$mut,side=1,at=1:n)
# par(op)

logger$info("Cumulative plots of regularization changes.")

cumuplot <- function(xs,percent=FALSE,add=FALSE,...) {
	xs <- sort(xs)
	max <- if (percent) 100 else 1
	# ys <- seq(0,max,length.out=length(xs))
	ys <- max * rank(xs)/length(xs)

	xlim <- list(...)[["xlim"]]
	if (!is.na(xlim[[1]])) {
		max.x <- xlim[[2]]
		xs <- c(xs,max.x)
		ys <- c(ys,max)
	}

	if (add) {
		lines(xs,ys,...)
	} else {
		plot(xs,ys,type="l",...)
	}
}

#How much change comes from the regularization?
geneNames <- c("UBE2I","SUMO1","TPK1","CALM1","NCS1")
geneColors <- c("orange","firebrick3","chartreuse3","steelblue3","violet")
names(geneColors) <- geneNames

html$subsection("Regularization score changes")
html$figure(function() {
	invisible(lapply(geneNames, function(geneName) {
		data <- read.csv(paste0(outdir,"imputed_regularized_",geneName,"_flipped_scores.csv"))
		changes <- with(data,na.omit(abs(screen.score-joint.score)))
		cumuplot(changes,percent=TRUE,
			add=(geneName!="UBE2I"),
			xlim=c(0,1),
			xlab=expression("regularization "~Delta~"score"),
			ylab=expression("% variants with "~Delta~"score < x"),
			col=geneColors[[geneName]]
		)
	}))
	legend("bottomright",geneNames,col=geneColors,lty=1)
},paste0(outdir,"regularization_change"),5,5)

html$shutdown()
logger$info("Done.")
