OUTDIR := workspace/$(shell date +"%Y%m%d-%H%M%S")/

all: scaleAndJoin

outdir:
	mkdir -p $(OUTDIR)

barseqTS: outdir
	Rscript bin/barseq-ts.R outdir=$(OUTDIR)

tileseq: outdir
	Rscript bin/tileseq.R outdir=$(OUTDIR) infile=input/raw_counts_SUMO1_tileseq.tsv geneName=SUMO1
	Rscript bin/tileseq.R outdir=$(OUTDIR) infile=input/raw_counts_UBE2I_tileseq.tsv geneName=UBE2I

scaleAndJoin: barseqTS tileseq
	Rscript bin/scaleAndJoin.R outdir=$(OUTDIR)
