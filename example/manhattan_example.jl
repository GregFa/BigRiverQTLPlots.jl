using BigRiverQTLPlots
using BulkLMM
using Random, Statistics
using Plots

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

@time single_results_perms = scan(
	pheno_y,
	geno_processed,
	kinship;
	permutation_test = true,
	nperms = 1000,
);

@time single_results = scan(
	pheno_y,
	geno_processed,
	kinship,
);


########
# Plot #
########

mycolor = ["#1f78b4", "#a6cee3", "#756bb1", "#bcbddc"];

plot_manhattan(
	single_results, gInfo,
	)
savefig(joinpath(@__DIR__, "..", "images", "manhattan_example.png"))

plot_manhattan(
	single_results_perms, gInfo, 
	significance = [0.10, 0.05],
	)
savefig(joinpath(@__DIR__, "..", "images", "manhattan_thrs_example.png"))
savefig(joinpath(@__DIR__, "..", "images", "manhattan_thrs_example.svg"))

plot_manhattan(
	single_results_perms, gInfo, 
	significance = [0.10, 0.05],
	manhattancolor = [mycolor[2] , mycolor[1]],
	)
savefig(joinpath(@__DIR__, "..", "images", "manhattan_thrs_example2.svg"))