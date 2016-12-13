options(stringsAsFactors=FALSE)

library("httr")
source("lib/libpdb.R")
source("lib/cliargs.R")
source("lib/liblogging.R")
source("lib/resultfile.R")
source("lib/trackdrawer.R")


pdb.acc <- getArg("accession",default="3uip")
main.chain <- getArg("chain",default="A")

outdir <- getArg("outdir",default="workspace/test/")


#Initialize logger
logger <- new.logger(paste0(outdir,"strucfeats-",pdb.acc,".log"))

#Set resultfile
html <- new.resultfile(paste0(outdir,"results.html"))
html$section(paste("Structural Feature extraction from",pdb.acc))


#delay timer
new.delay <- function(len=500) {
	lastMark <- NULL
	curr <- function() {
		as.numeric(Sys.time())*1000
	}
	elapsed <- function() {
		curr()-lastMark
	}
	mark <- function() {
		if (!is.null(lastMark) && elapsed() < len) {
			Sys.sleep((len-(elapsed()))/1000)
		}
		lastMark <<- curr()
	}
	list(elapsed=elapsed,mark=mark)
}


#query pdb for structure
query.pdb <- function(pdb.acc,pdb.file) {
	pdb.url <- paste0("https://files.rcsb.org/download/",pdb.acc,".pdb")
	htr <- GET(pdb.url)
	if (http_status(htr)$category == "Success") {
		writeBin(content(htr, "raw"), pdb.file)
	} else {
		logger$fatal("Unable to download PDB file")
		stop("Unable to download PDB file")
	}
}

#query getarea for solvent accessibility
query.getarea <- function(pdb.file) {
	getarea.url <- "http://curie.utmb.edu/cgi-bin/getarea.cgi"
	htr <- POST(getarea.url, body=list(
		water="1.4",
		gradient="n",
		name="popcode",
		email="jochen.weile@mail.utoronto.ca",
		Method="2",
		PDBfile=upload_file(pdb.file)
	))		
	if (http_status(htr)$category == "Success") {
		#parse getarea output
		con <- textConnection(content(htr,"text"))
		lines <- scan(con,what="character",sep="\n")
		close(con)
		lines <- lines[regexpr("----|POLAR|APOLAR|UNKNOW|Total|Number|\\*\\*",lines)<1]
		rawtable <- do.call(rbind,strsplit(lines," +"))
		return(data.frame(
			aa=rawtable[,2],
			pos=as.integer(rawtable[,3]),
			total.acc=as.numeric(rawtable[,4]),
			relative.acc=as.numeric(rawtable[,8])/100
		))
	} else {
		logger$fatal("Unable tor retrieve GETAREA results")
		stop("Unable to retrieve GETAREA results")
	}
}


#define filename for structure
pdb.file <- paste0(tempdir(),"/",pdb.acc,".pdb")

#read structure
logger$info("Querying PDB")
query.pdb(pdb.acc,pdb.file)
logger$info("Reading structure")
cplx.struc <- new.structure(pdb.file)
chains <- cplx.struc$get.chains()

logger$info("Splitting structure")
#Subdivide into chain pairs
subfiles <- cplx.struc$subcomplex.combos(main.chain)

#list of chain pairs
other.chains <- setdiff(chains,main.chain)
chain.sets <- c(main.chain,lapply(other.chains,c,main.chain))

delay <- new.delay(1000)

logger$info("Querying GETAREA")
#query getarea for all chain pairs
acctables <- mapply(function(cset,subfile) {
	#re-order chain ids to original file order
	cset <- intersect(chains,cset)
	ctags <- do.call(c,lapply(cset,function(x) rep(x,length(cplx.struc$get.sequence(x)))))
	delay$mark()#make sure to wait, as not to crowd the server
	cat("Querying subset",cset,"\n")
	acctable <- query.getarea(subfile)
	acctable[ctags==main.chain,]
},cset=chain.sets,subfile=subfiles,SIMPLIFY=FALSE)

logger$info("Computing burial")
#calculate burial based on accessibities
burial <- cbind(acctables[[1]],do.call(cbind,lapply(2:length(acctables),function(i) {
	before <- acctables[[1]][,"total.acc"]
	data.frame(
		abs.burial=before-acctables[[i]][,"total.acc"],
		rel.burial=(before-acctables[[i]][,"total.acc"])/(before+.000001)
	)
})))
colnames(burial)[seq(5,ncol(burial),2)] <- paste0("abs.burial.",other.chains)
colnames(burial)[seq(6,ncol(burial),2)] <- paste0("rel.burial.",other.chains)

query.stride <- function(pdb.file) {
	stride.url <- "http://webclu.bio.wzw.tum.de/cgi-bin/stride/stridecgi.py"
	htr <- POST(stride.url, body=list(
		inputfile=upload_file(pdb.file),
		action="compute"
	))
	if (http_status(htr)$category == "Success") {
		#parse getarea output
		con <- textConnection(content(htr,"text"))
		lines <- scan(con,what="character",sep="\n")
		close(con)
		lines <- lines[regexpr("ASG",lines)>0]
		rawtable <- do.call(rbind,strsplit(lines," +"))
		return(data.frame(
			aa=rawtable[,2],
			pos=as.integer(rawtable[,4]),
			secstruc=rawtable[,7]
		))
	} else {
		logger$fatal("Unable to retriev STRIDE results")
		stop("Unable to retrieve STRIDE results")
	}
}

logger$info("Querying STRIDE...")
secstruc <- query.stride(subfiles[[1]])

if (!all(secstruc$pos == burial$pos)) {
	logger$fatal("Inconsistent STRIDE results")
	stop("Inconsistent STRIDE results!")
}

logger$info("Writing results to file...")
outfile <- paste0(outdir,"strucfeats_",pdb.acc,".csv")
write.table(
	cbind(burial,ss=secstruc$secstruc),outfile,
	sep=",",row.names=FALSE,quote=FALSE
)
html$link.data(outfile)


burial <- cbind(burial,ss=secstruc$secstruc)
rownames(burial) <- burial$pos

allpos <- 1:max(burial$pos)
burial.all <- burial[as.character(allpos),]
rownames(burial.all) <- burial.all$pos <- allpos


html$figure(function() {
	td <- new.trackdrawer(l=nrow(burial.all))
	td$add.ss.track(burial.all$ss)
	td$add.track(burial.all$relative.acc,"Accessibility","steelblue3")
	for(oc in other.chains) {
		td$add.track(burial.all[,paste0("rel.burial.",oc)],oc,"orange",maxVal=1)
	}
	td$draw()
},paste0(outdir,"strucfeats_",pdb.acc),15,5)

logger$info("Done.")
