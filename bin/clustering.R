############
#Clustered genophenogram

options(stringsAsFactors=FALSE)
library("dendsort")
source("lib/cliargs.R")
source("lib/resultfile.R")
source("lib/liblogging.R")


outdir <- getArg("outdir",default="workspace/test/")


#Init logger
logger <- new.logger(paste0(outdir,"clustering.log"))

#Set resultfile
html <- new.resultfile(paste0(outdir,"results.html"))
html$section("Clustering genophenograms")

logger$info("Reading input")

proteins <- c("UBE2I","SUMO1","CALM1","TPK1","CBS")

scores <- do.call(rbind,lapply(proteins,function(prot) {
	scores <- read.csv(paste0(outdir,"imputed_regularized_",prot,"_flipped_scores.csv"))
	scores$pos <- with(scores,substr(mut,2,nchar(mut)-1))
	scores$maa <- with(scores,substr(mut,nchar(mut),nchar(mut)))
	scores$wtpos <- paste0(with(scores,substr(mut,1,nchar(mut)-1)),"-",prot)
	scores
}))

# scores <- read.csv(paste0(outdir,"imputed_regularized_UBE2I_scores.csv"))
# scores$pos <- with(scores,substr(mut,2,nchar(mut)-1))
# scores$maa <- with(scores,substr(mut,nchar(mut),nchar(mut)))
# scores$wtpos <- paste0(with(scores,substr(mut,1,nchar(mut)-1)),"-UBE2I")

# scores2 <- read.csv(paste0(outdir,"imputed_regularized_SUMO1_scores.csv"))
# scores2$pos <- with(scores2,substr(mut,2,nchar(mut)-1))
# scores2$maa <- with(scores2,substr(mut,nchar(mut),nchar(mut)))
# scores2$wtpos <- paste0(with(scores2,substr(mut,1,nchar(mut)-1)),"-SUMO1")

# scores <- rbind(scores,scores2)

logger$info("Aggregating")
aas <- unique(scores$maa)
score.by.pos <- with(scores,tapply(1:nrow(scores),wtpos,function(is) {
	if (!all(maa[is]==aas)) stop("wrong order!")
	joint.score[is]
}))

medians <- sapply(score.by.pos,median)

#Plot medians to decide on cutoff
op <- par(mfrow=c(length(proteins),1))
median.prots <- substr(names(medians),regexpr("-",names(medians))+1,nchar(names(medians)))
medians.by.prot <- tapply(medians,median.prots,c)
lapply(proteins,function(p) hist(
	medians.by.prot[[p]],
	main=p, xlab="positional medians",
	breaks=seq(-0.2,1,0.05),
	col="steelblue3",border=NA
))
par(op)

med.cutoff <- 0.5

lowscore.by.pos <- score.by.pos[medians < med.cutoff]
hiscore.by.pos <- score.by.pos[medians >= med.cutoff]

logger$info("Clustering")
pos.clust.lo <- dendsort(hclust(dist(do.call(rbind,lowscore.by.pos)),method="average"),type="average")
pos.clust.hi <- dendsort(hclust(dist(do.call(rbind,hiscore.by.pos)),method="average"),type="average")
# plot(pos.clust2)


geneColors <- c(UBE2I="gold",SUMO1="firebrick3",CALM1="steelblue3",TPK1="chartreuse3",CBS="purple")

logger$info("Drawing plot")

drawClusterGram <- function(pos.clust,score.by.pos,boxes=NULL,main="") {
	n <- length(pos.clust$order)
	posLabels <- pos.clust$labels[pos.clust$order]

	hIndex <- sapply(pos.clust$order,function(obs) which(apply(pos.clust$merge,1,function(pairing) -obs %in% pairing)))
	heights <- pos.clust$height[hIndex]

	layout(rbind(1,2))
	op <- par(mar=c(0,4,1,1)+.1,cex=0.7)
	plot(pos.clust,labels=FALSE,hang=0,main=main)
	wts <- substr(posLabels,1,1)
	text(1:n,min(heights)-0.1,wts,cex=0.7)
	segments(1:n,min(heights),1:n,heights,lty="dashed")
	par(las=1)
	plot(0,type="n",xlim=c(.5,n-.5),ylim=c(0,21),axes=FALSE,ylab="AA residue")
	colRamp <- colorRampPalette(c("royalblue3","white","firebrick3"))(11)
	colMatrix <- do.call(c,lapply(1:n, function(i) {
		pos <- posLabels[[i]]
		colIdx <- sapply(score.by.pos[[pos]],function(s) {
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
		colRamp[colIdx]
	}))
	x <- do.call(c,lapply(1:n,rep,20))
	y <- rep(20:1,n)
	rect(x-1,y-1,x,y,col=colMatrix,border=NA)

	sepLines <- c(7,9,12,15)
	segments(0,sepLines,n,sepLines,lty="dotted")
	arrows(-2,c(20,14.9,11.9,8.9),-2,c(15.1,12.1,9.1,7),length=0.01,angle=90,code=3)
	arrows(-4,c(20,11.9),-4,c(12.1,7),length=0.01,angle=90,code=3)
	text(-3,c(17.5,13.5,10.5,8),c(expression(psi),expression(phi),expression("+"),expression("-")))
	text(-5,c(16,9.5),c("hydrophobic","polar"),srt=90)

	axis(2,at=20:1-0.5,labels=aas)

	posGenes <- substr(posLabels,regexpr("-",posLabels)+1,nchar(posLabels))
	rect(1:n-1,20.5,1:n,21,col=geneColors[posGenes],border=NA)

	if (!is.null(boxes)) {
		rect(boxes[,1]-1,0,boxes[,2],20,lwd=2,lty="dashed")
	}

	par(op)
}

pdf("~/Dropbox/Documents/DMS/clustered_genophenogram.pdf",40,5)
	drawClusterGram(pos.clust.lo,lowscore.by.pos)
dev.off()

pdf("~/Dropbox/Documents/DMS/clustered_genophenogram_neutral.pdf",40,5)
	drawClusterGram(pos.clust.hi,hiscore.by.pos)
dev.off()

html$figure(function(){
	drawClusterGram(pos.clust2,lowscore.by.pos)
},paste0(outdir,"clustered_genophenogram"),20,5)


children <- function(clust) {
	ch <- as.list(rep(NA,nrow(clust$merge)))
	for (i in 1:nrow(clust$merge)) {
		left <- clust$merge[i,1]
		right <- clust$merge[i,2]
		leftch <- if (left < 0) -left else ch[[left]]
		rightch <- if (right < 0) -right else ch[[right]]
		ch[[i]] <- c(leftch,rightch)
	}
	ch
}

ideals <- function(chs, sbypos) {
	lapply(chs,function(ch) {
		rowMeans(do.call(cbind,sbypos[ch]))
	})
}

mses <- function(chs, sbypos, idls) {
	sapply(1:length(chs), function(i) {
		ch <- chs[[i]]
		idl <- idls[[i]]
		sum(sapply(sbypos[ch],function(sc) sum((sc-idl)^2)))/(20*length(ch))
	})
}

join.mse <- function(nodes, chs, mses) {
	if (any(nodes < 0)) {
		return(Inf)
	}
	ns <- sapply(chs[nodes],length)
	sum(mses[nodes]*20*ns)/(20*sum(ns))
}

children.lo <- children(pos.clust.lo)
ideals.lo <- ideals(children.lo, lowscore.by.pos)
mses.lo <- mses(children.lo, lowscore.by.pos, ideals.lo)

#greedy tree descend. split off clusters that allow greatest MSE decrease
# output: list of lists, with pairs "nodes" and "mse", where nodes is the list of cluster roots
descend.tree <- function(children.lo,mses.lo,pos.clust.lo) {
	root <- which.max(sapply(children.lo,length))
	descent <- list(
		list(nodes=root,mse=mses.lo[[root]]),
		list(nodes=pos.clust.lo$merge[root,],mse=join.mse(pos.clust.lo$merge[root,],children.lo,mses.lo))
	)
	for (ncuts in 3:200) {
		pnodes <- descent[[ncuts-1]]$nodes
		opts <- lapply(pnodes, function(cutnode) {
			kids <- pos.clust.lo$merge[cutnode,]
			nodes <- c(setdiff(pnodes,cutnode),kids)
			if (any(kids < 0)) {
				return(list(nodes=nodes,mse=Inf))
			} else {
				mse <- join.mse(nodes,children.lo,mses.lo)
				return(list(nodes=nodes,mse=mse))
			}
		})
		opt.mses <- sapply(opts,`[[`,"mse")
		if (all(is.infinite(opt.mses))) {
			cat("Singleton encountered at depth",ncuts)
			break;
		}
		best.opt <- which.min(opt.mses)
		descent[[ncuts]] <- opts[[best.opt]]
	}
	descent
}

descent <- descend.tree(children.lo,mses.lo,pos.clust.lo)

# op <- par(mfrow=c(2,1))
layout(rbind(1,2),heights=c(2,1))
plot(
	1:length(descent),
	sapply(descent,`[[`,"mse"),
	type="l", xlab="#clusters", ylab="MSE",
	ylim=c(0,0.13)
)
clSizes <- sapply(children.lo[descent[[72]]$nodes],length)
hist(clSizes,breaks=30,col="gray",border=NA,xlab="cluster size",main="72 clusters")
# par(op)

#########
# Re-draw clustergram with boxes for 72 clusters
boxes <- do.call(rbind,lapply(descent[[72]]$nodes,function(clroot) {
	range(sapply(children.lo[[clroot]],function(leaf) which(pos.clust.lo$order==leaf)))
}))
drawClusterGram(pos.clust.lo,lowscore.by.pos,boxes=boxes)


########
# Draw idealized archetypes for 72 clusters
clSizes.72 <- sapply(descent[[72]]$nodes,function(clroot) length(children.lo[[clroot]]) )
ideals.72 <- ideals.lo[descent[[72]]$nodes]
usefulness.72 <- sapply(ideals.72,function(x) sum(abs(x-0.5)))

clOrder <- order(usefulness.72,decreasing=TRUE)
clSizes.72 <- clSizes.72[clOrder]
clSizeCumu <- sapply(1:length(clSizes.72),function(i)sum(clSizes.72[1:i]))
ideals.72 <- ideals.72[clOrder]


op <- par(las=1,cex=0.7)
plot(0,type="n",xlim=c(0,sum(clSizes.72)),ylim=c(0,20),axes=FALSE,ylab="AA residue")
colRamp <- colorRampPalette(c("royalblue3","white","firebrick3"))(11)
rectcolors <- sapply(do.call(c,ideals.72), function(s) {
	colIdx <- {if (s <= 0) 1
	else if (s <= 0.2) 2
	else if (s <= 0.4) 3
	else if (s <= 0.6) 4
	else if (s <= .8) 5
	else if (s <= 1.2) 6
	else if (s <= 1.4) 7
	else if (s <= 1.6) 8
	else if (s <= 1.8) 9
	else if (s <= 2) 10
	else 11}
	colRamp[colIdx]
})
x1 <- do.call(c,lapply(c(0,clSizeCumu[-length(clSizeCumu)]),rep,20))
x2 <- do.call(c,lapply(clSizeCumu,rep,20))
y <- rep(20:1,length(ideals.72))
rect(x1,y-1,x2,y,col=rectcolors,border=NA)

sepLines <- c(7,9,12,15)
segments(0,sepLines,sum(clSizes.72),sepLines,lty="dotted")
arrows(-2,c(20,14.9,11.9,8.9),-2,c(15.1,12.1,9.1,7),length=0.01,angle=90,code=3)
arrows(-4,c(20,11.9),-4,c(12.1,7),length=0.01,angle=90,code=3)
text(-3,c(17.5,13.5,10.5,8),c(expression(psi),expression(phi),expression("+"),expression("-")))
text(-5,c(16,9.5),c("hydrophobic","polar"),srt=90)

axis(2,at=20:1-0.5,labels=aas)
par(op)




######
# Draw clustergrams for each AA 
#
pdf("~/Dropbox/Documents/DMS/clustersByAA.pdf",11,8.5)
invisible(lapply(aas, function(aa) {
	score.by.pos.k <- score.by.pos[medians < med.cutoff & sapply(score.by.pos,function(v)v[[which(aas==aa)]]>0.7)]
	if (length(score.by.pos.k) > 0) {
		pos.clust.k <- dendsort(hclust(dist(do.call(rbind,score.by.pos.k)),method="average"),type="average")
		drawClusterGram(pos.clust.k,score.by.pos.k,main=aa)
	}
}))
dev.off()



alphabet <- read.csv("~/alphaPoly.csv")
lettercolors <- c(
	A="steelblue1",V="steelblue2",L="steelblue3",I="steelblue",M="steelblue4",
	F="cyan1",Y="cyan2",W="cyan3",R="firebrick1",H="firebrick2",K="firebrick3",
	D="maroon1",E="maroon2",S="chartreuse1",T="chartreuse2",N="chartreuse3",Q="chartreuse4",
	G="orange",C="salmon",P="gold"
)
drawLetter <- function(letter, x, y, w, h) {
	polygon(
		x=w*with(alphabet,x[pathGroup==paste0(letter,".1")])+x,
		y=h*with(alphabet,y[pathGroup==paste0(letter,".1")])+y,
		col=lettercolors[[letter]],border=NA
	)
	if (any(alphabet$pathGroup == paste0(letter,".2"))) {
		polygon(
			x=w*with(alphabet,x[pathGroup==paste0(letter,".2")])+x,
			y=h*with(alphabet,y[pathGroup==paste0(letter,".2")])+y,
			col="white",border=NA
		)
	}
	if (any(alphabet$pathGroup == paste0(letter,".3"))) {
		polygon(
			x=w*with(alphabet,x[pathGroup==paste0(letter,".3")])+x,
			y=h*with(alphabet,y[pathGroup==paste0(letter,".3")])+y,
			col="white",border=NA
		)
	}
}

pdf("~/Dropbox/Documents/DMS/logo.pdf",40,5)
plot(0,type="n",xlim=c(0,72),ylim=c(0,1),axes=FALSE,xlab="",ylab="")
j <- 0
invisible(lapply(ideals.72,function(vals) {
	vals <- sapply(vals,function(x) if(x < 0) 0 else if (x < 0.5) x/3 else x)
	relVals <- vals/sum(vals)
	cumRelVals <- 1-sapply(1:20,function(i)sum(relVals[1:i]))
	for (i in 1:20) {
		if (relVals[[i]] > 0.02) {
			drawLetter(letter=aas[[i]],x=j,y=,cumRelVals[[i]],w=1,h=relVals[[i]])
		}
	}
	j <<- j+1
}))
dev.off()


#permute
permutation.data <- do.call(rbind,replicate(100,{
	perm.scores <- lapply(lowscore.by.pos,sample,20)
	#cluster
	pos.clust.perm <- hclust(dist(do.call(rbind,perm.scores)),method="average")
	#leaves for each node
	children.perm <- children(pos.clust.perm)
	nleaves.perm <- sapply(children.perm,length)
	#idealized profile for each node
	ideals.perm <- ideals(children.perm, perm.scores)
	#MSEs for each node
	mses.perm <- mses(children.perm, perm.scores, ideals.perm)
	data.frame(nleaves=nleaves.perm,mse=mses.perm)
},simplify=FALSE))
permutation.distros <- with(permutation.data,tapply(mse,nleaves,c))


#Plot some example distributions
op <- par(mfrow=c(5,1))
for (i in 2:6) {
	distro <- permutation.distros[[as.character(i)]]
	hist(distro,xlim=c(0,0.15),
		main=paste("Cluster size =",i),
		freq=FALSE,xlab="MSE",
		breaks=seq(0,0.2,0.002),
		col="gray",border=NA
	)
	# lines(seq(0,0.15,0.001),dnorm(seq(0,0.15,0.001),mean(distro),sd(distro)),col=2)
}
par(op)


mses.72 <- mses.lo[descent[[72]]$nodes]
pvals.72 <- sapply(1:72,function(i) {
	csize <- clSizes.72[[i]]
	cmse <- mses.72[[i]]
	distro <- permutation.distros[[as.character(csize)]]
	pnorm(cmse,mean(distro),sd(distro))
})


# lowscore.by.pos
# pos.clust.lo
# children.lo
# ideals.lo
# mses.lo

#filter down to top-level sets (i.e. remove subclusters already contained in other clusters)
topset <- function(is,leaves) {
	hasUnique <- sapply(is,function(i) {
		all(sapply(setdiff(is,i),function(j) {
			#any unique leaves in i?
			length(setdiff(leaves[[i]],leaves[[j]])) > 0
		}))
	})
	return(is[hasUnique])
}


clSizesAll <- sapply(children.lo,length)
pvals.all <- sapply(1:length(children.lo), function(i) {
	csize <- clSizesAll[[i]]
	cmse <- mses.lo[[i]]
	if (as.character(csize) %in% names(permutation.distros)) {
		distro <- permutation.distros[[as.character(csize)]]
		# pnorm(cmse,mean(distro),sd(distro))
		sum(distro <= cmse)/length(distro)
	} else {
		NA
	}
})

candidates <- which(clSizesAll < length(lowscore.by.pos)/3 & clSizesAll > 2)
qvals <- p.adjust(pvals.all[candidates])
qvals.all <- 1
qvals.all[candidates] <- qvals

# pvals.all.adj <- p.adjust(pvals.all,method="fdr")

# signif.clust <- which(qvals < 0.1 & clSizesAll < length(lowscore.by.pos)/3)
signif.clust <- which(qvals.all < 0.1)
signif.clust.top <- topset(signif.clust,children.lo)
# for (i in 1:length(signif.clust)) {
# 	for (j in i:length(signif.clust)) {
# 		c.i <- children.lo[[signif.clust[[i]]]]
# 		c.j <- children.lo[[signif.clust[[j]]]]
# 		cat(i,",",j,"->",length(intersect(c.i,c.j)),"\n")
# 	}
# }

#########
# Re-draw clustergram with boxes for significant clusters
boxes <- do.call(rbind,lapply(signif.clust,function(clroot) {
	range(sapply(children.lo[[clroot]],function(leaf) which(pos.clust.lo$order==leaf)))
}))
drawClusterGram(pos.clust.lo,lowscore.by.pos,boxes=boxes)




######
# Draw clustergrams for each WT 
#
wt.dendros <- lapply(aas, function(aa) {
	wt.aas <- substr(names(score.by.pos),1,1)
	score.by.pos.k <- score.by.pos[medians < med.cutoff & wt.aas==aa]
	if (length(score.by.pos.k) > 0) {
		pos.clust.k <- dendsort(hclust(dist(do.call(rbind,score.by.pos.k)),method="average"),type="average")
		leaves.k <- children(pos.clust.k)
		ideals.k <- ideals(leaves.k, score.by.pos.k)
		mses.k <- mses(leaves.k, score.by.pos.k, ideals.k)
		csizes.k <- sapply(leaves.k,length)
		pvals.k <- mapply(function(mse,csize) {
			distro <- permutation.distros[[as.character(csize)]]
			sum(distro <= mse)/length(distro)
		},mse=mses.k,csize=csizes.k)
		# drawClusterGram(pos.clust.k,score.by.pos.k,main=aa)
		return(list(scores=score.by.pos.k,dendro=pos.clust.k,
			leaves=leaves.k,ideals=ideals.k,mses=mses.k,pvals=pvals.k
		))
	} else {
		return(NULL)
	}
})
names(wt.dendros) <- aas

pdf("~/Dropbox/Documents/DMS/clustersByWT.pdf",11,8.5)
invisible(lapply(aas, function(aa) {
	wt.dendro <- wt.dendros[[aa]]
	if (length(wt.dendro) > 0) {
		with(wt.dendro,{
			qvals <- p.adjust(pvals)
			signifs <- which(qvals < 0.05)
			boxes <- NULL
			if (length(signifs) > 0) {
				topsignifs <- topset(signifs,leaves)
				boxes <- do.call(rbind,lapply(signifs,function(signif) {
					range(sapply(leaves[[signif]],function(leaf) which(dendro$order==leaf)))
				}))
			}
			drawClusterGram(dendro,scores,main=aa,boxes=boxes)
		})
	}
}))
dev.off()




html$shutdown()


logger$info("Done!")
