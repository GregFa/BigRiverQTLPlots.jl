
using BigRiverQTLPlots
using BulkLMM
using Random, Statistics
using Plots
using Helium

##############
# BXD spleen #
##############

########
# Data #
########
bulklmmdir = dirname(pathof(BulkLMM));

gmap_file = joinpath(bulklmmdir, "..", "data", "bxdData", "gmap.csv");
gInfo = BulkLMM.CSV.read(gmap_file, BulkLMM.DataFrames.DataFrame);

phenocovar_file = joinpath(bulklmmdir, "..", "data", "bxdData", "phenocovar.csv");
pInfo = BulkLMM.CSV.read(phenocovar_file, BulkLMM.DataFrames.DataFrame);

pheno_file = joinpath(bulklmmdir, "..", "data", "bxdData", "spleen-pheno-nomissing.csv");
pheno = BulkLMM.DelimitedFiles.readdlm(pheno_file, ',', header = false);
pheno_processed = pheno[2:end, 2:(end-1)] .* 1.0; # exclude the header, the first (transcript ID)and the last columns (sex)

geno_file = joinpath(bulklmmdir, "..", "data", "bxdData", "spleen-bxd-genoprob.csv")
geno = BulkLMM.DelimitedFiles.readdlm(geno_file, ',', header = false);
geno_processed = geno[2:end, 1:2:end] .* 1.0;

#################
# Preprocessing #
#################
traitID = 1112;
pheno_y = reshape(pheno_processed[:, traitID], :, 1);

###########
# Kinship #
###########
kinship = calcKinship(geno_processed);

########
# Scan #
########

single_results_perms = scan(
	pheno_y,
	geno_processed,
	kinship;
	permutation_test = true,
	nperms = 1000,
);

single_results = scan(
	pheno_y,
	geno_processed,
	kinship,
);

thrs = BigRiverQTLPlots.perms_thresholds(single_results_perms.L_perms, [0.10, 0.05]);
# Helium.writehe(reshape(thrs, :, 1), joinpath(@__DIR__, "..", "test", "thresholds.he"))
########
# Plot #
########
plot_QTL(single_results.lod, gInfo)
savefig(joinpath(@__DIR__, "..", "images", "QTL_example.png"))

plot_QTL(single_results_perms, gInfo, significance = [0.10, 0.05])
savefig(joinpath(@__DIR__, "..", "images", "QTL_thrs_example.png"))

