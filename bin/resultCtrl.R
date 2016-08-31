
source("lib/resultfile.R")
source("lib/cliargs.R")

outdir <- getArg("outdir",default="workspace/test/")
cmd <- getArg("cmd",default="setup")

html <- new.resultfile(paste0(outdir,"results.html"))

if (cmd=="setup") {
	html$header("POPCode pipeline")
} else {
	html$closing()
}

html$shutdown()
