############
#Clustered genophenogram

options(stringsAsFactors=FALSE)
library("dendsort")
source("lib/cliargs.R")
source("lib/resultfile.R")
source("lib/liblogging.R")


outdir <- getArg("outdir",default="workspace/test/")


#Init logger
logger <- new.logger(paste0(outdir,"geneticInteractions.log"))

#Set resultfile
html <- new.resultfile(paste0(outdir,"results.html"))
html$section("Clustering genophenograms")

logger$info("Reading input")

scores <- read.csv(paste0(outdir,"imputed_regularized_UBE2I_scores.csv"))
scores$pos <- with(scores,substr(mut,2,nchar(mut)-1))
scores$maa <- with(scores,substr(mut,nchar(mut),nchar(mut)))
scores$wtpos <- paste0(with(scores,substr(mut,1,nchar(mut)-1)),"-UBE2I")

scores2 <- read.csv(paste0(outdir,"imputed_regularized_SUMO1_scores.csv"))
scores2$pos <- with(scores2,substr(mut,2,nchar(mut)-1))
scores2$maa <- with(scores2,substr(mut,nchar(mut),nchar(mut)))
scores2$wtpos <- paste0(with(scores2,substr(mut,1,nchar(mut)-1)),"-SUMO1")

scores <- rbind(scores,scores2)

logger$info("Clustering")
aas <- unique(scores$maa)
score.by.pos <- with(scores,tapply(1:nrow(scores),wtpos,function(is) {
	if (!all(maa[is]==aas)) stop("wrong order!")
	joint.score[is]
}))
pos.clust2 <- dendsort(hclust(dist(do.call(rbind,score.by.pos)),method="average"),type="average")
# plot(pos.clust2)
n <- length(pos.clust2$order)

logger$info("Drawing plot")

# pdf("clustered_genophenogram.pdf",14,5)
html$figure(function(){
	layout(rbind(1,2))
	op <- par(mar=c(0,4,1,1)+.1)
	plot(pos.clust2,cex=.7)
	plot(0,type="n",xlim=c(0,n+1),ylim=c(0,20),axes=FALSE)
	colRamp <- colorRampPalette(c("royalblue3","white","firebrick3"))(11)
	invisible(lapply(1:n, function(i) {
		pos <- pos.clust2$labels[pos.clust2$order[[i]]]
		colIdx <- sapply(score.by.pos[[pos.clust2$labels[pos.clust2$order[[i]]]]],function(s) {
			if (s <= 0) 1
			else if (s <= 0.2) 2
			else if (s <= 0.4) 3
			else if (s <= 0.6) 4
			else if (s <= .8) 5
			else if (s <= 1.2) 6
			else if (s <= 1.4) 7
			else if (s <= 1.6) 8
			else if (s <= 1.8) 9
			else if (s <= 2) 10
			else 11
		})
		cols <- colRamp[colIdx]
		rect(i-1,19:0,i,20:1,col=cols,border=NA)
	}))
},paste0(outdir,"clustered_genophenogram"),23,5)
# dev.off()

html$shutdown()


logger$info("Done!")
