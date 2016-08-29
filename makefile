OUTDIR := workspace/$(shell date +"%Y%m%d-%H%M%S")/

all: evaluateSpotting compensatory pickSpottingClones

#Create an output directory for this pipeline run
outdir:
	mkdir -p $(OUTDIR)

#Analyze the BarSEQ time series data from raw barcode counts
barseqTS: outdir
	Rscript bin/barseq-ts.R outdir=$(OUTDIR)

#Analyze the TileSEQ complementation data from raw mutation counts
tileseq: outdir
	Rscript bin/tileseq.R outdir=$(OUTDIR) infile=input/raw_counts_SUMO1_tileseq.tsv geneName=SUMO1
	Rscript bin/tileseq.R outdir=$(OUTDIR) infile=input/raw_counts_UBE2I_tileseq.tsv geneName=UBE2I

#Re-scale the TileSEQ data to match the BarSEQ scale and join the datasets
scaleAndJoin: barseqTS tileseq
	Rscript bin/scaleAndJoin.R outdir=$(OUTDIR)

#Build a feature table and then use RandomForest prediction to impute missing data
# and regularize poorly measured data points
impute: scaleAndJoin
	Rscript bin/impute.R outdir=$(OUTDIR) infile=$(OUTDIR)compl_joint_results_UBE2I.csv \
		geneName=UBE2I
	Rscript bin/impute.R outdir=$(OUTDIR) infile=$(OUTDIR)compl_tileSEQ_results_SUMO1_transformed.csv \
		geneName=SUMO1

#Use double-mutant information to detect intragenic epistasis
geneticInteractions: scaleAndJoin
	Rscript bin/geneticInteractions.R outdir=$(OUTDIR) infile=$(OUTDIR)compl_joint_results_UBE2I.csv

#Use crystal structure data to calculate inter-C_alpha distances of all amino acid pairs
distanceMatrix: outdir
	Rscript bin/distanceMatrix.R outdir=$(OUTDIR)

#Use intragenic epistasis and C_alpha distances to find candidate compensatory relationships
compensatory: geneticInteractions distanceMatrix
	Rscript bin/compensatoryMut.R outdir=$(OUTDIR) infile=$(OUTDIR)genetic_interactions.csv \
		distfile=$(OUTDIR)distanceMatrix_UBE2I.csv

#Pick a set of clones to be tested in a manual spotting assay
pickSpottingClones: scaleAndJoin
	Rscript bin/pickSpottingClones.R outdir=$(OUTDIR)

#Compare the results of the manual spotting assay to the results of the
# previous pipeline steps (barseq, tileseq, impute, etc)
evaluateSpotting: impute
	Rscript bin/evaluateSpotting.R outdir=$(OUTDIR)
