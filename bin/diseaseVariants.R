#######################################
# Assess disease variants #
#######################################
options(stringsAsFactors=FALSE)

source("lib/resultfile.R")
source("lib/liblogging.R")
source("lib/cliargs.R")
source("lib/libroc.R")
source("lib/libyogitools.R")

library("beeswarm")


outdir <- getArg("outdir",default="workspace/test/")
# outdir <- getArg("outdir",default="workspace/20170131-175115/")

#Initialize logger
logger <- new.logger(paste0(outdir,"diseaseVariants.log"))

#Set resultfile
html <- new.resultfile(paste0(outdir,"results.html"))
html$section("Disease Variant Analysis")


geneNames <- c("CALM1","CALM2","CALM3","NCS1","TPK1","UBE2I","SUMO1")
proteinNames <- c("CALM1","CALM1","CALM1","NCS1","TPK1","UBE2I","SUMO1")


# extract.groups <- function(x, re) {
# 	matches <- regexpr(re,x,perl=TRUE)
# 	start <- attr(matches,"capture.start")
# 	end <- start + attr(matches,"capture.length") - 1
# 	do.call(cbind,lapply(1:ncol(start), function(i) {
# 		sapply(1:nrow(start),function(j){
# 			if (start[j,i] > -1) substr(x[[j]],start[j,i],end[j,i]) else NA
# 		})
# 	}))
# }

read.gnomad <- function(infile,geneName) {
	aa.letters <- c(
		Ala="A",Val="V",Ile="I",Leu="L",Met="M",Phe="F",Tyr="Y",Trp="W",
		Arg="R",His="H",Lys="K",Asp="D",Glu="E",Ser="S",Thr="T",Asn="N",Gln="Q",
		Cys="C",Gly="G",Pro="P"
	) 
	rawtable <- read.csv(infile)
	filtertable <- rawtable[
		with(rawtable,Annotation=="missense" & regexpr("†",Consequence)<0),
		c("Protein.Consequence","Allele.Frequency","Number.of.Homozygotes")
	]
	with(filtertable,data.frame(
		gene=geneName,
		db="gnomad",
		mut=sapply(Protein.Consequence,function(s) {
			from <- aa.letters[substr(s,3,5)]
			to <- aa.letters[substr(s,nchar(s)-2,nchar(s))]
			pos <- substr(s,6,nchar(s)-3)
			paste0(from,pos,to)
		}),
		freq=Allele.Frequency,
		homozygotes=Number.of.Homozygotes
	))
}

gnomad <- do.call(rbind,lapply(
	geneNames,
	function(gn)read.gnomad(paste0("res/alleles/gnomad_",gn,".csv"),gn)
))


read.cosmic <- function(infile,geneName) {
	rawtable <- read.csv(infile)
	filtertable <- rawtable[regexpr("Missense",rawtable$Type)>0,]
	with(filtertable,data.frame(
		gene=geneName,
		db="cosmic",
		mut=sub("p.","",AA.Mutation),
		freq=Count
	))
}

geneNames <- c("CALM1","CALM2","CALM3","NCS1","TPK1","UBE2I","SUMO1")

cosmic <- do.call(rbind,lapply(
	geneNames,
	function(gn)read.cosmic(paste0("res/alleles/cosmic_",gn,".csv"),gn)
))


read.clinvar <- function(infile,geneName) {
	aa.letters <- c(
		Ala="A",Val="V",Ile="I",Leu="L",Met="M",Phe="F",Tyr="Y",Trp="W",
		Arg="R",His="H",Lys="K",Asp="D",Glu="E",Ser="S",Thr="T",Asn="N",Gln="Q",
		Cys="C",Gly="G",Pro="P"
	) 
	rawtable <- read.delim(infile)
	filtertable <- rawtable[rawtable$Gene.s.==geneName,]
	if (nrow(filtertable)==0) return(NULL)
	# extract.groups(filtertable$Name,"\\(p\\.(\\w{3})(\\d+)(\\w{3})\\)")
	filtertable <- filtertable[!is.na(extract.groups(filtertable$Name,"\\(p\\.(\\w{3})(\\d+)(\\w{3})\\)")[,1]),]
	if (nrow(filtertable)==0) return(NULL)
	muts <- apply(extract.groups(filtertable$Name,"\\(p\\.(\\w{3})(\\d+)(\\w{3})\\)"),1,function(x) {
		paste0(aa.letters[x[[1]]],x[[2]],aa.letters[x[[3]]])
	})
	anno <- sub("\\(.+\\)","",filtertable$Clinical.significance..Last.reviewed)
	data.frame(
		gene=geneName,
		db="clinvar",
		mut=muts,
		annotation=anno
	)
}

geneNames <- c("CALM1","CALM2","CALM3","NCS1","TPK1","UBE2I","SUMO1")

clinvar <- do.call(rbind,lapply(
	geneNames,
	function(gn) {
		# cat(gn,"\t")
		read.clinvar(paste0("res/alleles/clinvar_",gn,".txt"),gn)
	}
))

####
#Kick out any variants that do not reference the correct amino acids
####
validateAA <- function(muts,geneNames) {
	parseFASTA <- function(filename) scan(filename,what="character",sep="\n")[[2]]
	seqs <- c(
		CALM1=parseFASTA("res/CALM1_aa.fa"),
		CALM2=parseFASTA("res/CALM1_aa.fa"),
		CALM3=parseFASTA("res/CALM1_aa.fa"),
		NCS1=parseFASTA("res/NCS1_aa.fa"),
		TPK1=parseFASTA("res/TPK1_aa.fa"),
		UBE2I=parseFASTA("res/UBE2I_aa.fa"),
		SUMO1=parseFASTA("res/SUMO1_aa.fa")
	)
	mapply(function(mut,gene){
		pos <- as.integer(substr(mut,2,nchar(mut)-1))
		aa <- substr(mut,1,1)
		substr(seqs[[gene]],pos,pos)==aa
	},mut=muts,gene=geneNames)
}

gnomad <- gnomad[with(gnomad,validateAA(mut,gene)),]
cosmic <- cosmic[with(cosmic,validateAA(mut,gene)),]
clinvar <- clinvar[with(clinvar,validateAA(mut,gene)),]

# read.ltri <- function(filename) {
# 	data <- read.csv(filename)
# 	with(data,data.frame(row.names=mut,score=joint.score,sd=joint.sd))
# }

scores <- lapply(proteinNames,function(prot) {
	filename <- paste0(outdir,"imputed_regularized_",prot,"_flipped_scores.csv")
	data <- read.csv(filename)
	with(data,data.frame(row.names=mut,score=joint.score,se=joint.se))
})
names(scores) <- geneNames


fetchScores <- function(muts,geneNames) {
	
	do.call(rbind,mapply(function(mut,gene){
		# aa <- substr(mut,nchar(mut),nchar(mut))
		# pos <- as.integer(substr(mut,2,nchar(mut)-1))
		# gpgs[[gene]][aa,pos]
		scores[[gene]][mut,]
	},mut=muts,gene=geneNames,SIMPLIFY=FALSE))

}

gnomad <- cbind(gnomad,with(gnomad,fetchScores(mut,gene)))
cosmic <- cbind(cosmic,with(cosmic,fetchScores(mut,gene)))
clinvar <- cbind(clinvar,with(clinvar,fetchScores(mut,gene)))

outfile <- paste0(outdir,"gnomad_vars.csv")
write.table(gnomad,outfile,sep=",",row.names=FALSE)
html$link.data(outfile)

outfile <- paste0(outdir,"cosmic_vars.csv")
write.table(cosmic,outfile,sep=",",row.names=FALSE)
html$link.data(outfile)

outfile <- paste0(outdir,"clinvar_vars.csv")
write.table(clinvar,outfile,sep=",",row.names=FALSE)
html$link.data(outfile)



kgenomes <- data.frame(
	gene="TPK1",db="1000Genomes",
	mut=c("A145E","C192R","E57K","G223W","L17W","R33H","R51H","V226D","Y31S"),
	mendel="heterozygous",score=1
)

html$subsection("Variant distributions")
html$figure(function() {

	op <- par(mfrow=c(3,2)) 

	breaks <- seq(-4,1,0.1)
	hist(scores[["CALM1"]][,"score"],breaks=breaks,xlim=c(-0.5,1),col="gray",border=NA,main="CALM1/2/3",xlab="score")
	abline(v=with(gnomad,score[regexpr("CALM",gene)>0]),col="chartreuse4")
	abline(v=with(clinvar,score[regexpr("CALM",gene)>0]),col="firebrick3")
	legend("topleft",
		c("rare polymorphism","common polymorphism","pathogenic allele","benign allele",
			"singleton somatic tumor mutation","somatic tumor mutation (n>1)"),
		col=c("chartreuse4","chartreuse1","firebrick3","blue","gold","orange"),
		lty="solid"
	)

	hist(scores[["TPK1"]][,"score"],breaks=breaks,xlim=c(-0.5,1),col="gray",border=NA,main="TPK1 diploid",xlab="score")
	abline(v=with(kgenomes,score[gene=="TPK1"]),col="chartreuse4")
	# abline(v=with(gnomad,score[gene=="TPK1"]),col=sapply(with(gnomad,freq[gene=="TPK1"])<1e-3,ifelse,"chartreuse4","chartreuse1"))
	abline(v=with(clinvar,score[gene=="TPK1"]),col=sapply(with(clinvar,annotation[gene=="TPK1"])=="Benign",ifelse,"blue","firebrick3"))

	hist(scores[["TPK1"]][,"score"],breaks=breaks,xlim=c(-0.5,1),col="gray",border=NA,main="TPK1 haploid",xlab="score")
	abline(v=with(gnomad,score[gene=="TPK1"]),col=sapply(with(gnomad,freq[gene=="TPK1"])<1e-3,ifelse,"chartreuse4","chartreuse1"))
	abline(v=with(clinvar,score[gene=="TPK1"]),col=sapply(with(clinvar,annotation[gene=="TPK1"])=="Benign",ifelse,"blue","firebrick3"))

	hist(scores[["NCS1"]][,"score"],breaks=breaks,xlim=c(-0.5,1),col="gray",border=NA,main="NCS1",xlab="score")
	abline(v=with(gnomad,score[gene=="NCS1"]),col="chartreuse4")
	abline(v=with(cosmic,score[gene=="NCS1"]),col=sapply(with(cosmic,freq[gene=="NCS1"])>1,ifelse,"orange","gold"))
	
	hist(scores[["UBE2I"]][,"score"],breaks=breaks,xlim=c(-0.5,1),col="gray",border=NA,main="UBE2I",xlab="score")
	abline(v=with(gnomad,score[gene=="UBE2I"]),col="chartreuse4")
	abline(v=with(cosmic,score[gene=="UBE2I"]),col=sapply(with(cosmic,freq[gene=="UBE2I"])>1,ifelse,"orange","gold"))
	
	hist(scores[["SUMO1"]][,"score"],breaks=breaks,xlim=c(-0.5,1),col="gray",border=NA,main="SUMO1",xlab="score")
	abline(v=with(gnomad,score[gene=="SUMO1"]),col="chartreuse4")
	abline(v=with(cosmic,score[gene=="SUMO1"]),col=sapply(with(cosmic,freq[gene=="SUMO1"])>1,ifelse,"orange","gold"))
	
	par(op)

},paste0(outdir,"diseaseVariants"),14,7)

##########
#converts a number in scientific notation into an expression for plotting
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

html$subsection("Somatic variants vs Polymorphisms")
html$figure(function() {
	swarm <- list(
		GnomAD=with(gnomad,score[gene %in% c("UBE2I","SUMO1","NCS1")]),
		COSMIC=with(cosmic,score[gene %in% c("UBE2I","SUMO1","NCS1")])
	)
	beeswarm(swarm,pch=20,col=c("chartreuse4","orange"),
		ylim=c(-.5,1.5),labels=c("polymorphisms","somatic tumor\nmutations"),
		ylab="score"
	)
	abline(h=0:1,col=c("firebrick3","chartreuse3"))
	bxplot(swarm,add=TRUE,col="gray40")
	pval <- with(swarm,wilcox.test(GnomAD,COSMIC,alternative="greater"))$p.value
	if (pval < 0.05) {
		lines(c(1,1,2,2),c(1.1,1.2,1.2,1.1))
		text(1.5,1.3,scinotExpr(pval))
	}
},paste0(outdir,"somaticVsPoly"),5,5)


html$subsection("Somatic variants vs Polymorphisms (w/o NCS1)")
html$figure(function() {
	swarm <- list(
		GnomAD=with(gnomad,score[gene %in% c("UBE2I","SUMO1")]),
		COSMIC=with(cosmic,score[gene %in% c("UBE2I","SUMO1")])
	)
	beeswarm(swarm,pch=20,col=c("chartreuse4","orange"),
		ylim=c(-.5,1.5),labels=c("polymorphisms","somatic tumor\nmutations"),
		ylab="score"
	)
	abline(h=0:1,col=c("firebrick3","chartreuse3"))
	bxplot(swarm,add=TRUE,col="gray40")
	pval <- with(swarm,wilcox.test(GnomAD,COSMIC,alternative="greater"))$p.value
	if (pval < 0.05) {
		lines(c(1,1,2,2),c(1.1,1.2,1.2,1.1))
		text(1.5,1.3,scinotExpr(pval))
	}
},paste0(outdir,"somaticVsPoly_noNCS1"),5,5)

html$subsection("Precision-Recall curves")
html$figure(function() {

	op <- par(mfrow=c(2,1))
	truth <- c(
		rep(FALSE,sum(regexpr("CALM",gnomad$gene)>0)),
		rep(TRUE,sum(regexpr("CALM",clinvar$gene)>0))
	)
	pred <- c(
		1-with(gnomad,score[regexpr("CALM",gene)>0]),
		1-with(clinvar,score[regexpr("CALM",gene)>0])
	)
	roc.obj <- roc.data(truth,pred)
	draw.prc(roc.obj,main="CALM1/2/3")
	auc.atlas <- auprc(roc.obj)

	provean <- read.delim("res/alleles/CALM_provean.txt")
	roc.obj2 <- roc.data(truth,10-provean$SCORE)
	draw.prc(roc.obj2,col="blue",add=TRUE)
	auc.provean <- auprc(roc.obj2)

	polyphen <- read.csv("res/alleles/CALM_polyphen.csv")
	roc.obj3 <- roc.data(truth,polyphen$pph2_prob)
	draw.prc(roc.obj3,col="chartreuse3",add=TRUE)
	auc.polyphen <- auprc(roc.obj3)

	legend("bottomleft",
		c(
			sprintf("DMS Atlas AUC = %.02f",auc.atlas),
			sprintf("PROVEAN AUC = %.02f",auc.provean),
			sprintf("Polyphen-2 AUC = %.02f",auc.polyphen)
		),
		lwd=2,col=c("firebrick3","blue","chartreuse3")
	)

	truth <- c(
		rep(FALSE,sum(kgenomes$gene=="TPK1")),
		clinvar[clinvar$gene=="TPK1","annotation"]!="Benign"
	)
	pred <- c(
		1-with(kgenomes,score[gene=="TPK1"]),
		1-with(clinvar,score[gene=="TPK1"])
	)
	roc.obj <- roc.data(truth,pred)
	draw.prc(roc.obj,main="TPK1")

	provean <- read.delim("res/alleles/TPK1_provean.tsv")
	provean.pred <- c(rep(0,sum(kgenomes$gene=="TPK1")),10-provean$SCORE)
	roc.obj2 <- roc.data(truth,provean.pred)
	draw.prc(roc.obj2,col="blue",add=TRUE)

	polyphen <- read.delim("res/alleles/TPK1_polyphen.tsv")
	pp.pred <- c(rep(0,sum(kgenomes$gene=="TPK1")),polyphen$pph2_prob)
	roc.obj3 <- roc.data(truth,pp.pred)
	draw.prc(roc.obj3,col="chartreuse3",add=TRUE)

	par(op)

},paste0(outdir,"vus_prc"),5,10)


html$shutdown()
logger$info("Done.")
