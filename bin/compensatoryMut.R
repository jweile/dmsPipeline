
options(stringsAsFactors=FALSE)

# source("lib/libyogitools.R")
source("lib/liblogging.R")
source("lib/cliargs.R")

library(beeswarm)

outdir <- getArg("outdir",default="workspace/test/")
infile <- getArg("infile",default="workspace/test/genetic_interactions.csv")

#Init logger
logger <- new.logger(paste0(outdir,"compensatoryMut.log"))

##############
# LOAD INPUT #
##############
logger$info("Reading input")

dist <- read.csv("res/ube2i_distances.csv")
colnames(dist) <- rownames(dist)

gis <- read.csv(infile)

mutpos <- do.call(rbind,lapply(strsplit(gis$dm,","),function(ms) as.integer(substr(ms,2,nchar(ms)-1))))

##################
# INTEGRATE DATA #
##################
logger$info("Integrating distance information")

gis$dist <- apply(mutpos,1,function(x) {
	p1 <- as.character(x[[1]])
	p2 <- as.character(x[[2]])
	if (!all(c(p1,p2) %in% colnames(dist))) {
		return(NA)
	}
	dist[p1,p2]
})

logger$info("Saving integrated data table")
write.table(gis,paste0(outdir,"gi+dist.csv"),sep=",",row.names=FALSE)

###############
# FILTER DATA #
###############
logger$info("Filter for compensatory candidates")
compens <- gis[which(with(gis,epsilon > 0.2 & dist < 10)),]


################
# DRAW FIGURES #
################
logger$info("Drawing figure")

pdf(paste0(outdir,"compensatory.pdf"),5,5)
op <- par(mar=c(5,5,4,1)+.1,las=1)
with(gis,plot(epsilon,dist,
	xlab=expression(epsilon[i][j]),
	ylab=expression(scriptstyle(bgroup("||",list(C[alpha]^(i),C[alpha]^(j)),"||")[2])),
	pch=16,col=rgb(79,148,205,100,maxColorValue=255),
	main="Genetic interaction vs Residue distance",axes=FALSE,
	ylim=c(0,50)
))
axis(1)
axis(2,at=seq(0,50,10),labels=paste(seq(0,50,10),"Å"))
# rect(0,10,3,-1,border="firebrick3",density=6,col="firebrick3",lty="dashed")
rect(0.2,10,3,-1,col=rgb(200,20,20,50,maxColorValue=255),border="firebrick3",lty="dotted")
text(1.3,4.5,"compensatory zone",col="firebrick3")
with(compens,text(epsilon,dist,sub(",","-",dm),cex=.5,srt=30,pos=4))
par(op)
invisible(dev.off())

detail <- function(i) {
	vals <- c(WT=1,s1=compens[i,4],s2=compens[i,6],d=compens[i,2])
	sds <- c(WT=NA,s1=compens[i,5],s2=compens[i,7],d=compens[i,3])
	ses <- sapply(sds,function(x)sqrt(x)/6)
	labs <- c("WT",strsplit(compens[i,1],",")[[1]],"DM")
	xs <- barplot(vals,names.arg=labs,
		col="lightgoldenrod3",border="gray",
		ylab="Complementation fitness"
	)
	arrows(xs[,1],vals-ses/2,xs[,1],vals+ses/2,angle=90,length=0.05,code=3)
}

pdf(paste0(outdir,"compensatory_detail.pdf"),4,8)
layout(rbind(1,2))
detail(1)
detail(2)
invisible(dev.off())

logger$info("Done.")
