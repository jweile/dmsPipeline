####################################
# Pick representative clones for   #
# manual spotting assay            #
####################################

# source("lib/libyogitools.R")
source("lib/liblogging.R")
source("lib/cliargs.R")
# library("hash")

options(stringsAsFactors=FALSE)

#get output directory
outdir <- getArg("outdir",default="workspace/test/")

#Initialize logger
logger <- new.logger(paste0(outdir,"pickSpottingClones.log"))

##############
# LOAD INPUT #
##############
logger$info("Reading input")

joint.data <- read.csv(paste0(outdir,"compl_joint_results_UBE2I.csv"))
joint.muts <- strsplit(joint.data$mut,",")
joint.nmut <- sapply(joint.muts,length)
rownames(joint.data) <- joint.data$mut

bc.clones <- read.csv("res/clones.csv")
rownames(bc.clones) <- bc.clones$id

#Update clone table with deletion information
muts <- strsplit(bc.clones$aa.calls,",")
nmut <- sapply(muts,length)

delpos <- do.call(rbind,lapply(strsplit(bc.clones$deletions,"-"),function(x) if (length(x)==1&&is.na(x)) c(NA,NA) else as.numeric(x)))
foo <- delpos[,1] < 82 & delpos[,2] > 545
foo[is.na(foo)] <- FALSE
null.clones <- bc.clones$id[foo]
foo <- delpos[,1] < delpos[,2]
foo[is.na(foo)] <- FALSE
longdel.clones <- bc.clones$id[foo]

bc.clones$final.call <- bc.clones$aa.calls
bc.clones[longdel.clones,"final.call"] <- "longdel"
bc.clones[null.clones,"final.call"]<- "null"

#Pick single mutant clones with high genotype confidence
good.singles <- with(bc.clones,bc.clones[!(final.call %in% c("WT","longdel","null")) & freq > 0.6 & nmut == 1,])
#and make a list of the corresponding genotypes
bc.single.muts <- unique(good.singles$final.call)

logger$info("Finding imputation control set")

#Find the set of single AA mutations that is not represented in the screen results
unrepresented <- setdiff(bc.single.muts,joint.data$mut)
#and find the matching clones
unrep.clones <- good.singles[good.singles$final.call %in% unrepresented,]

logger$info("Picking random score spectrum clones")

#Find the subset of genotypes with good screen data quality (=high-fidelity set)
hifi <- na.omit(with(joint.data,joint.data[mut %in% bc.single.muts & sd < .2 & joint.nmut == 1 & regexpr("_",mut) < 1,]))

#Define 64 bins linearly across the score spectrum in the high-fidelity set
bins <- seq(min(hifi$score),max(hifi$score),length.out=64)

#For each bin, pick one high-quality mutation at random that falls into its interval.
spectrum.muts <- na.omit(sapply(1:(length(bins)-1),function(i) {
	muts <- with(hifi,mut[score >= bins[[i]] & score < bins[[i+1]]])
	if (length(muts) > 0) sample(muts,1) else NA
}))
# length(spectrum.muts)

#Pick clones that correspond to the spectrum genotypes
spectrum.clones <- good.singles[good.singles$final.call %in% spectrum.muts,]
#Remove any duplicates
spectrum.clones <- spectrum.clones[!duplicated(spectrum.clones$final.call),]
#And look up out the associatioed screen scores
spectrum.scores <- joint.data[spectrum.clones$final.call,"score"]

logger$info("Writing results to file.")

#Combine the results in one output table
out <- data.frame(
	category=c(rep("imputation",length(unrepresented)),rep("spectrum",length(spectrum.muts))),
	id=c(unrep.clones$id,spectrum.clones$id),
	plate=c(unrep.clones$dest.plate,spectrum.clones$dest.plate),
	well=c(unrep.clones$dest.well,spectrum.clones$dest.well),
	mut=c(unrep.clones$final.call,spectrum.clones$final.call),
	score=c(rep(NA,length(unrepresented)),spectrum.scores)
)

write.table(out,paste0(outdir,"picked_spottingAssay_clones.csv"),sep=",",row.names=FALSE)

logger$info("Done.")
