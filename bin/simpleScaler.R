
source("lib/resultfile.R")
source("lib/liblogging.R")
source("lib/cliargs.R")

options(stringsAsFactors=FALSE)


#get output directory
outdir <- getArg("outdir",default="workspace/test/")
#get input file
infile <- getArg("infile",default="workspace/test/compl_tileSEQ_results_UBE2I.csv")
#get gene name
geneName <- getArg("geneName",default="UBE2I")

#Initialize logger
logger <- new.logger(paste0(outdir,"simpleScaler-",geneName,".log"))

#Set resultfile
html <- new.resultfile(paste0(outdir,"results.html"))
html$section(paste(geneName,"Simple TileSEQ scaler"))

##############
# LOAD INPUT #
##############
logger$info("Reading input")

data <- read.csv(infile)

transf <- with(data,{
	stop.is <- which(substr(mut,nchar(mut),nchar(mut))=="_")
	syn.is <- which(substr(mut,1,1)==substr(mut,nchar(mut),nchar(mut)))
	# miss.is <- setdiff(1:nrow(data),c(stop.is,syn.is))

	med.stop <- median(mean.lphi[stop.is],na.rm=TRUE)
	med.syn <- median(mean.lphi[syn.is],na.rm=TRUE)
	denom <- med.syn - med.stop
	sd <- bsd/denom
	cbind(
		score=(mean.lphi - med.stop)/denom,
		sd=sd,
		df=df,
		se=sd/sqrt(df)
	)
})

outfile <- paste0(outdir,"compl_scaled_results_",geneName,".csv")
write.table(data.frame(mut=data$mut,transf),outfile,sep=",",row.names=FALSE)
html$link.data(outfile)

html$shutdown()

logger$info("Done.")
