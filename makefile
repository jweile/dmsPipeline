OUTDIR := workspace/$(shell date +"%Y%m%d-%H%M%S")/

all: impute compensatory pickSpottingClones

outdir:
	mkdir -p $(OUTDIR)

barseqTS: outdir
	Rscript bin/barseq-ts.R outdir=$(OUTDIR)

tileseq: outdir
	Rscript bin/tileseq.R outdir=$(OUTDIR) infile=input/raw_counts_SUMO1_tileseq.tsv geneName=SUMO1
	Rscript bin/tileseq.R outdir=$(OUTDIR) infile=input/raw_counts_UBE2I_tileseq.tsv geneName=UBE2I

scaleAndJoin: barseqTS tileseq
	Rscript bin/scaleAndJoin.R outdir=$(OUTDIR)

impute: scaleAndJoin
	Rscript bin/impute.R outdir=$(OUTDIR) infile=$(OUTDIR)compl_joint_results_UBE2I.csv \
		geneName=UBE2I
	Rscript bin/impute.R outdir=$(OUTDIR) infile=$(OUTDIR)compl_tileSEQ_results_SUMO1_transformed.csv \
		geneName=SUMO1

geneticInteractions: scaleAndJoin
	Rscript bin/geneticInteractions.R outdir=$(OUTDIR) infile=$(OUTDIR)compl_joint_results_UBE2I.csv

distanceMatrix: outdir
	Rscript bin/distanceMatrix.R outdir=$(OUTDIR)

compensatory: geneticInteractions distanceMatrix
	Rscript bin/compensatoryMut.R outdir=$(OUTDIR) infile=$(OUTDIR)genetic_interactions.csv \
		distfile=$(OUTDIR)distanceMatrix_UBE2I.csv

pickSpottingClones: scaleAndJoin
	Rscript bin/pickSpottingClones.R outdir=$(OUTDIR)
