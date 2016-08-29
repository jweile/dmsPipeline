
source("lib/libpdb.R")
source("lib/liblogging.R")
source("lib/cliargs.R")

#get output directory
outdir <- getArg("outdir",default="workspace/test/")

#Initialize logger
logger <- new.logger(paste0(outdir,"distanceMatrix.log"))

logger$info("Reading input")

#Load the UBE2I structure from PDB entry 3UIP
ube2i <- new.structure("res/3uip.pdb")

#Get the list of all AA positions in the UBE2I structure
positions <- as.integer(names(ube2i$get.sequence("A")))

#Define a function for Euclidian Distance
eucl.dist <- function(x,y) sqrt(sum((x-y)^2))

logger$info("Calculating distance matrix")

#Calculate the matrix of Euclidian distances between all alpha carbon pairs
pb <- txtProgressBar(max=length(positions),style=3)
pr <- 0
dmat <- do.call(rbind,lapply(positions, function(p.i) {
	setTxtProgressBar(pb,value=pr<<-pr+1)
	coord.i <- unlist(ube2i$get.atom("A",p.i,"CA")[,c("x","y","z"),drop=TRUE])
	sapply(positions, function(p.j) {
		coord.j <- unlist(ube2i$get.atom("A",p.j,"CA")[,c("x","y","z"),drop=TRUE])
		eucl.dist(coord.i,coord.j)
	})
}))
close(pb)

dimnames(dmat) <- list(positions,positions)


#Draw heatmap figure

logger$info("Drawing figure")

pdf(paste0(outdir,"distanceMatrix.pdf"))
layout(cbind(1,2),widths=c(8,1.5))
#draw distance map
op <- par(mar=c(5,4,4,0)+.1)
fire <- colorRampPalette(c(c("white","yellow","red","black")))(25)
image(dmat,axes=FALSE,col=fire,xlab="AA position",ylab="AA position")
labels <- 1:15*10
at <- which(positions %in% labels)/length(positions)
axis(1,at=at,labels=labels)
axis(2,at=at,labels=labels)
#draw legend
par(mar=c(5,0,4,4)+.1)
plot(0,type="n",xlim=0:1,ylim=c(0,25),axes=FALSE,xlab="",ylab="")
rect(0,0:24,1,1:25,col=fire,border=NA)
axis(4,at=seq(0.5,24.5,length.out=25),
	labels=round(seq(min(dmat),max(dmat),length.out=25))
)
mtext(expression(C[alpha]~"distance (Å)"),side=4,line=3)
par(op)
invisible(dev.off())

#Output

logger$info("Writing output to file.")

write.table(dmat,paste0(outdir,"distanceMatrix_UBE2I.csv"),sep=",")

logger$info("Done.")
