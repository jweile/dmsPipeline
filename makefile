OUTDIR := workspace/$(shell date +"%Y%m%d-%H%M%S")/

all: evaluateSpotting colorizeStructure findInterfaces compensatory accCons yeastResidues finalize

#Create an output directory for this pipeline run
outdir:
	mkdir -p $(OUTDIR)
	Rscript bin/resultCtrl.R outdir=$(OUTDIR) cmd=setup

#Analyze the BarSEQ time series data from raw barcode counts
barseqTS: outdir
	Rscript bin/barseq-ts.R outdir=$(OUTDIR)

#Analyze the TileSEQ complementation data from raw mutation counts
tileseq: outdir
	Rscript bin/tileseq.R outdir=$(OUTDIR) infile=input/raw_counts_SUMO1_tileseq.tsv geneName=SUMO1
	Rscript bin/tileseq.R outdir=$(OUTDIR) infile=input/raw_counts_UBE2I_tileseq.tsv geneName=UBE2I

#Analyze the BarSEQ time series data from raw barcode counts
barseqY2H: outdir
	Rscript bin/barseq-y2h.R outdir=$(OUTDIR)

#Re-scale the TileSEQ data to match the BarSEQ scale and join the datasets
scaleAndJoin: barseqTS tileseq
	Rscript bin/scaleAndJoin.R outdir=$(OUTDIR)

#Build a feature table and then use RandomForest prediction to impute missing data
# and regularize poorly measured data points
impute: scaleAndJoin
	Rscript bin/impute.R outdir=$(OUTDIR) infile=$(OUTDIR)compl_joint_results_UBE2I.csv \
		geneName=UBE2I ctrlSet=res/UBE2I_imputationControl.txt
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
evaluateSpotting: impute pickSpottingClones
	Rscript bin/evaluateSpotting.R outdir=$(OUTDIR)

#Plot the relationship between mutant fitness and surface accessibility,
# interfacialness and conservation
accCons: impute
	Rscript bin/accCons.R outdir=$(OUTDIR)

colorizeStructure: impute
	Rscript bin/colorizeStructure.R outdir=$(OUTDIR) \
		infile=$(OUTDIR)imputed_regularized_UBE2I_scores.csv geneName=UBE2I
	Rscript bin/colorizeStructure.R outdir=$(OUTDIR) \
		infile=$(OUTDIR)imputed_regularized_SUMO1_scores.csv geneName=SUMO1

findInterfaces: barseqY2H impute
	Rscript bin/findInterfaces.R outdir=$(OUTDIR)

yeastResidues: impute
	Rscript bin/yeastResidues.R outdir=$(OUTDIR)

#Adds closing tags to the result HTML
finalize: outdir
	Rscript bin/resultCtrl.R outdir=$(OUTDIR) cmd=finalize
