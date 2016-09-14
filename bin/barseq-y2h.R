####################################
# Process Raw mutation counts from #
# Y2H barseq data                  #
####################################

source("lib/resultfile.R")
source("lib/liblogging.R")
source("lib/libyogitools.R")
source("lib/cliargs.R")

library("hash")
options(stringsAsFactors=FALSE)


quantileNorm <- as.logical(getArg("quantileNorm",default=FALSE))
if (is.na(quantileNorm)) {
	quantileNorm <- FALSE
}

#get output directory
outdir <- getArg("outdir",default="workspace/test/")

#Initialize logger
logger <- new.logger(paste0(outdir,"barseq-y2h.log"))

#Set resultfile
html <- new.resultfile(paste0(outdir,"results.html"))
html$section("UBE2I Y2H BarSEQ analysis")


logger$info("Reading input")

sample.table <- read.delim("res/muxtags_y2h_samples.tsv")
data.mat <- read.csv("input/raw_counts_Y2H.csv")
# clones <- rownames(data.mat)

data.mat2 <- data.mat[-1,paste("X",sample.table$sample,sep="")]

data.rel <- apply(data.mat2,2,function(x)x/sum(x))


quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
   
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
   
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}
#apply quantile normalization to replicate groups
if (quantileNorm) {
	logger$info("Performing quantile normalization")
	invisible(lapply(with(sample.table,tapply(sample,paste(interactor,condition,sep="_"),c)), function(is) {
		data.rel[,is] <<- quantile_normalisation(data.rel[,is])
	}))
}

logger$info("Plotting barcode frequency distributions")

#BOXPLOT
html$subsection("Barcode frequency distributions")
html$figure(function(){
	ord <- with(sample.table,order(interactor,condition,replicate))
	#more colors: "peachpuff","chartreuse","aquamarine","darkgoldenrod","lightsteelblue"
	plotcols <- do.call(c,(lapply(do.call(c,lapply(c("firebrick","darkolivegreen","steelblue"),function(x) paste(x,c(1,3,4),sep=""))),rep,3)))
	boxplot(data.rel[,ord],ylim=c(0,.0001),col=plotcols,ylab="rel.freq.",axes=FALSE)
	axis(2)
	invisible(with(sample.table[ord,],tapply(1:nrow(sample.table),interactor,function(is) {
		mtext(interactor[is[[1]]],side=1,line=2,at=mean(is))
		tapply(is,condition[is],function(js) {
			mtext(condition[js[[1]]],side=1,line=1,at=mean(js))
		})
	})))
},paste0(outdir,"y2h_freqdist"),10,5)


############################
###MEASURE JACKPOT EFFECT###
############################
logger$info("Plotting jackpot effect diagnostic")

html$subsection("Jackpot effect diagnostic")
html$figure(function(){
	ord <- with(sample.table,order(interactor,condition,replicate))	
	plotcols <- do.call(c,(lapply(do.call(c,lapply(c("firebrick","darkolivegreen","steelblue"),function(x) paste(x,c(1,3,4),sep=""))),rep,3)))
	colpos <- barplot(apply(data.rel[,ord],2,max),col=plotcols,
		ylab="max(rel.freq.)",axes=FALSE,names.arg=NA,
		main="Jackpot effect"
	)
	axis(2)
	invisible(with(sample.table[ord,],tapply(1:nrow(sample.table),interactor,function(is) {
		mtext(interactor[is[[1]]],side=1,line=2,at=mean(colpos[is,]))
		tapply(is,condition[is],function(js) {
			mtext(condition[js[[1]]],side=1,line=1,at=mean(colpos[js,]))
		})
	})))
},paste0(outdir,"y2h_jackpot"),10,5)


logger$info("Loading clone library")

#Load clones
clone.table <- read.csv("res/clones_y2h.csv",stringsAsFactors=FALSE)
clone.idx <- hash(clone.table$id,1:nrow(clone.table))
muts <- strsplit(clone.table$aa.calls,",")


delpos <- do.call(rbind,lapply(strsplit(clone.table$deletions,"-"),function(x) if (length(x)==1&&is.na(x)) c(NA,NA) else as.numeric(x)))
foo <- delpos[,1] < 82 & delpos[,2] > 545
foo[is.na(foo)] <- FALSE
null.clones <- clone.table$id[foo]
foo <- delpos[,1] < delpos[,2]
foo[is.na(foo)] <- FALSE
longdel.clones <- clone.table$id[foo]
foo <- delpos[,1] > delpos[,2]
foo[is.na(foo)] <- FALSE
dupl.clones <- clone.table$id[foo]


logger$info("Calculating fitness values for each strain")

interactors <- unique(sample.table$interactor)
all.data <- NULL

# pdf("techrep_y2hV4.pdf",12,4)
html$subsection("Technical replication")
html$figure(function(){
	op <- par(mfrow=c(length(interactors),3))
	all.data <<- do.call(cbind,lapply(interactors,function(ia) {

		permissive.cols <- with(sample.table,sample[interactor==ia & condition=="+HIS"])
		selective.cols <- with(sample.table,sample[interactor==ia & condition=="-HIS"])
		medBC <- apply(data.mat2[,permissive.cols],1,median)
		well.measured <- medBC > 100

		logfc <- log(apply(data.rel[,selective.cols],1,mean) / apply(data.rel[,permissive.cols],1,mean))

		logfc.rep <- log(data.rel[,selective.cols] / data.rel[,permissive.cols])
		sds <- apply(logfc.rep,1,sd)

		logfc.rep.plot <- do.call(rbind,lapply(1:nrow(logfc.rep),function(i) {
			if (any(is.na(logfc.rep[i,])) || any(is.infinite(logfc.rep[i,]))) NULL else logfc.rep[i,]
		}))

		invisible(apply(combn(3,2),2,function(is) {
			plot(logfc.rep.plot[,is],
				xlim=c(-7,5),ylim=c(-7,5),pch=".",col="steelblue3",
				xlab=paste("lfc rep",is[[1]]),ylab=paste("lfc rep",is[[2]]),
				main=paste(ia,is[[1]],"vs",is[[2]])
			)
			text(0,4,paste("R =",signif(cor(na.omit(logfc.rep.plot[,is]))[1,2],3)))
		}))
		# hist(logfc,breaks=50,col="darkolivegreen3",border="gray40")

		ztab <- data.frame(log.medBC=log10(medBC),lfc=logfc)
		deadrows <- which(medBC==0 | is.na(sds) | sds==0)
		z <- lm(log10(sds[-deadrows])~., ztab[-deadrows,])
		model.sd <- 10^predict(z,ztab)
		#Baldi&Long
		bnl <- function(pseudo.n,n,model.sd,empiric.sd) {
			sqrt((pseudo.n * model.sd^2 + (n - 1) * empiric.sd^2)/(pseudo.n + n - 2))
		}
		bsd <- bnl(6,3,model.sd,sds)


		out <- data.frame(lfc=logfc,medBC=medBC,esd=sds,bsd=bsd)
		colnames(out) <- paste(ia,c("lfc","medBC","esd","bsd"),sep=".")
		out

	}))
	par(op)
},paste0(outdir,"y2h_techrep"),6,8)
# dev.off()


html$subsection("Permissive condition barcode distribution")
html$figure(function(){
	with(all.data,hist(
		log10(SATB1.medBC),breaks=50,
		col="gray",border=NA,axes=FALSE,
		xlab="median +HIS barcode count",
		main=""
	))
	abline(v=2,col="red")
	axis(2)
	axis(1,at=0:6,labels=10^(0:6))
},paste0(outdir,"y2h_wmcounts"),10,5)



# add genotype to data
all.data$mut <- clone.table[values(clone.idx,rownames(data.rel)),"aa.calls"]
all.data[longdel.clones,"mut"] <- "longdel"
all.data[null.clones,"mut"] <- "null"
all.data[dupl.clones,"mut"] <- "longdup"
all.data$single <- regexpr("null|longdel|,",all.data$mut) < 0

#Filter out uncertain genotypes and indels
logger$info("Filtering out bad genotypes")
good.clones <- with(clone.table,id[sapply(freq,function(x) is.na(x) || x > .6) & call != "longdel" & call != "longdup"])
good.clones <- good.clones[good.clones %in% rownames(all.data)]
good.data <- all.data[good.clones,]


#BLOSUM and Polyphen
# library("Biostrings")
# data(BLOSUM62)

# good.data$min.blosum <- sapply(1:nrow(good.data), function(i) {
# 	# if (good.data$single[[i]]) {
# 		ms <- good.data$mut[[i]]
# 		if (!(ms %in% c("null","longdel"))) {
# 			min(sapply(ms, function(m) {
# 				from <- substr(m,1,1)
# 				to <- substr(m,nchar(m),nchar(m))
# 				BLOSUM62[from,to]
# 			}))
# 		} else NA
# 	# } else NA
# })

# pp.in <- read.delim("res/polyphen_ube2i.txt",comment.char="#",stringsAsFactors=FALSE,header=FALSE)
# polyphen <- with(pp.in,data.frame(row.names=paste(V8,V7,V9,sep=""),pp=V12))

# single.data <- good.data[good.data$single,]
# single.data$pp <- polyphen[single.data$mut,"pp"]


#Extract null and WT controls
null.clones <- rownames(good.data)[with(good.data,which(regexpr("NULL",rownames(good.data))>0 & SATB1.medBC > 100))]
wt.clones <- rownames(good.data)[with(good.data,which(regexpr("WT",rownames(good.data))>0 & SATB1.medBC > 100))]

logger$info("Plotting score distributions")

html$subsection("Well-measured score distribution with controls")
html$figure(function(){
	layout(rbind(1,2,3,4))
	invisible(lapply(interactors,function(ia) {

		null.lfc <- na.omit(good.data[null.clones,paste(ia,"lfc",sep=".")])
		null.lfc <- null.lfc[is.finite(null.lfc)]
		wt.lfc <- na.omit(good.data[wt.clones,paste(ia,"lfc",sep=".")])
		wt.lfc <- wt.lfc[is.finite(wt.lfc)]
		wm <- good.data[,paste(ia,"medBC",sep=".")] > 100
		hist(good.data[wm,paste(ia,"lfc",sep=".")],breaks=50,col="gray",border=NA,main=paste(ia,"48h"),xlab="LFC")
		abline(v=null.lfc,col="firebrick3")
		abline(v=wt.lfc,col="darkolivegreen3")
		abline(v=median(wt.lfc,na.rm=TRUE),col="darkolivegreen3",lwd=5,lty="dashed")
		abline(v=median(null.lfc),col="firebrick3",lwd=5,lty="dashed")
		good.data[,paste(ia,"i.score",sep=".")] <<- (good.data[,paste(ia,"lfc",sep=".")] - median(null.lfc))/(median(wt.lfc)-median(null.lfc))
		#TODO: Also re-scale stdev
	}))
},paste0(outdir,"y2h_wmscores"))



logger$info("Plotting biological replication")

html$subsection("Biological replication")
html$figure(function(){
	op <- par(mfrow=c(2,2))
	invisible(lapply(interactors,function(ia) {
		
		sac.data <- good.data[!(good.data$mut %in% c("null","wt","WT","longdel","longdup")) & good.data[,paste(ia,"medBC",sep=".")]>100,]

		#lfc pairs for biological replicates
		br.pairs <- fin2(do.call(rbind,tapply(sac.data[,paste(ia,"i.score",sep=".")],sac.data$mut,function(x) {
			if (length(x) > 1) t(combn(x,2)) else NULL
		},simplify=FALSE)))
		plot(br.pairs,xlim=c(-1,2),ylim=c(-1,2),pch=20,main=ia,col="steelblue3",
			xlab="score biol. repl. 1",ylab="score biol. repl. 2")
		text(0.5,1.5,paste("R =",signif(cor(br.pairs)[1,2],3)))
	}))
	par(op)
},paste0(outdir,"y2h_biorep"))


logger$info("Writing scores to file")

outfile <- paste0(outdir,"y2h_scores.csv")
write.table(good.data,outfile,sep=",")
html$link.data(outfile)

html$shutdown()

logger$info("Done.")

