#######################################
# Assess somatic and natural variants #
#######################################

source("lib/resultfile.R")
# source("lib/libyogitools.R")
source("lib/liblogging.R")
# source("lib/topoScatter.R")
source("lib/cliargs.R")
# library("hash")
library("beeswarm")

tstat <- function (n1,n2,m1,m2,sd1,sd2) {
  tt <- -(m1 - m2)/sqrt((((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2)/(n1 + n2 - 2)) * ((n1 + n2)/(n1 * n2)))
  dft <- n1 + n2 - 2
  # t^2 follows F distribution
  p <- 1 - pf(tt^2, 1, dft)
  c(t=tt,p=p)
}

options(stringsAsFactors=FALSE)

#get output directory
outdir <- getArg("outdir",default="workspace/test/")
#get gene name
# geneName <- getArg("geneName",default="UBE2I")

#Initialize logger
logger <- new.logger(paste0(outdir,"somaticVnatural-",geneName,".log"))

#Set resultfile
html <- new.resultfile(paste0(outdir,"results.html"))
html$section(paste(geneName,"Somatic vs natural variants"))


soma.nat <- function(geneName) {

	logger$info("Reading input")
	screen <- read.csv(paste0(outdir,"imputed_regularized_",geneName,"_scores.csv"))
	rownames(screen) <- screen$mut

	logger$info("Calculating health scores")

	#Convert complementation score into "health" score by treating
	# significantly hyperactive mutants as deleterious
	wt.sd <- 0.402 #(taken from WT control stdev in barseq)
	health.q <- p.adjust(sapply(1:nrow(screen),function(i) {
		s <- screen$joint.score[[i]]
		if (s > 1) {
			sd <- screen$joint.sd[[i]]
			tstat(10,10,1,s,wt.sd,sd)[["p"]]
		} else NA
	}),method="fdr")
	screen$healthy <- sapply(1:nrow(screen),function(i){
		if (!is.na(health.q[[i]]) && health.q[[i]] < 0.05){
			1/screen$joint.score[[i]]^2
		} else
		screen$joint.score[[i]]
	})

	screen.missense <- screen[sapply(screen$mut,function(m)substr(m,1,1)!=substr(m,nchar(m),nchar(m))),]

	logger$info("Reading natural variant data")

	aa.letters <- c(
		Ala="A",Val="V",Ile="I",Leu="L",Met="M",Phe="F",Tyr="Y",Trp="W",
		Arg="R",His="H",Lys="K",Asp="D",Glu="E",Ser="S",Thr="T",Asn="N",Gln="Q",
		Cys="C",Gly="G",Pro="P"
	)

	natural <- read.csv(paste0("res/",geneName,"_gnomad.csv"))
	natural <- natural[natural$Annotation=="missense",]
	natural.muts <- data.frame(
		mut=sapply(natural$Protein.Consequence,function(s) {
			from <- substr(s,3,5)
			to <- substr(s,nchar(s)-2,nchar(s))
			pos <- substr(s,6,nchar(s)-3)
			paste0(aa.letters[from],pos,aa.letters[to])
		}),
		occurrence=sapply(natural$Allele.Frequency,function(f) {
			if (f < 1e-6) {
				"very.rare"
			} else if (f < 1/1000) {
				"rare"
			} else "common"
		})
	)
	natural.muts <- natural.muts[which(natural.muts$mut %in% screen$mut),]
	natural.muts$fitness <- screen[natural.muts$mut,"healthy"]


	logger$info("Reading somatic variant data")

	somatic <- read.csv(paste0("res/",geneName,"_cosmic.csv"))
	somatic <- somatic[somatic$Type=="Substitution - Missense",]
	somatic.muts <- data.frame(
		mut=sub("p.","",somatic$AA.Mutation),
		occurrence=sapply(somatic$Count,function(s) if (s < 2) "singleton" else "multiple")
	)
	somatic.muts <- somatic.muts[which(somatic.muts$mut %in% screen$mut),]
	somatic.muts$fitness <- screen[somatic.muts$mut,"healthy"]

	vars <- rbind(
		cbind(class=rep("natural",nrow(natural.muts)),natural.muts),
		cbind(class=rep("somatic",nrow(somatic.muts)),somatic.muts)
	)
	rownames(vars) <- NULL
	list(vars=vars,all=screen.missense$healthy)
}

ubdata <- soma.nat("UBE2I")
sudata <- soma.nat("SUMO1")

allvars <- rbind(ubdata$vars,sudata$vars)
allfits <- c(ubdata$all,sudata$all)


draw.vars <- function(allvars,allfits) {
	layout(cbind(1,2),widths=c(2,1))
	hist(allfits,breaks=100,col="gray",border=NA,
		xlab="fitness",main="",xlim=c(-0.5,1.6),prob=TRUE
	)
	with(allvars,abline(v=fitness[class=="natural"],
		col="gray40",lwd="2",lty="dashed"
	))
	with(allvars,abline(v=fitness[class=="somatic"],
		col=sapply(occurrence[class=="somatic"],
			function(o) if(o=="singleton")"orange" else "firebrick3"
		),
		lwd="2",lty="dashed"
	))
	legend("topright",c("natural","somatic singleton","somatic multiple"),
		col=c("gray40","orange","firebrick3"),lwd=2,lty="dashed",bg="white"
	)

	# fitness <- c(somatic.muts$fitness,natural.muts$fitness)
	# variants <- c(rep("somatic",nrow(somatic.muts)),rep("natural",nrow(natural.muts)))
	beecol <- sapply(allvars$occurrence,
		function(s)
		if (s == "rare") "gray40"
		else if (s == "singleton") "orange"
		else if (s == "multiple") "firebrick3"
		else "black"
	)
	# beecol <- c(sapply(somatic.muts$singleton,ifelse,"orange","firebrick3"),rep("gray40",nrow(natural.muts)))
	maxy <- max(allvars$fitness)

	with(allvars,{
		bxplot(fitness ~ class,col="gray",ylab="fitness",ylim=c(-0.5,maxy+0.5))
		beeswarm(fitness ~ class,pwcol=beecol,pch=20,cex=1.5,add=TRUE)
		abline(h=0:1,col=c("firebrick3","chartreuse3"))

		wilcox <- wilcox.test(fitness[class=="somatic"], fitness[class=="natural"],alternative="less")
		if (wilcox$p.value < 0.05) {
			lines(c(1,1,2,2),c(max(fitness[class=="natural"])+.1,maxy+0.2,maxy+0.2,max(fitness[class=="somatic"])+.1))
			text(1.5,maxy+0.3,"*",cex=1.5)
		}
	})
}



logger$info("Plotting variant distribution")

html$subsection("UBE2I variants")
html$figure(function(){
	draw.vars(ubdata$vars,ubdata$all)
},paste0(outdir,"somaticVnatural_",geneName),10,5)

html$subsection("SUMO1 variants")
html$figure(function(){
	draw.vars(sudata$vars,sudata$all)
},paste0(outdir,"somaticVnatural_",geneName),10,5)

html$subsection("UBE2I & SUMO1 variants")
html$figure(function(){
	draw.vars(allvars,allfits)
},paste0(outdir,"somaticVnatural_",geneName),10,5)


html$shutdown()
logger$info("Done!")

