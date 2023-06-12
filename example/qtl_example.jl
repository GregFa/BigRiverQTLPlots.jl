
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

# idx_geno = findall(occursin.(gInfo.Chr, "1 2 3 4 5"));
# gInfo_subset = gInfo[idx_geno, :];

phenocovar_file = joinpath(bulklmmdir, "..", "data", "bxdData", "phenocovar.csv");
pInfo = BulkLMM.CSV.read(phenocovar_file, BulkLMM.DataFrames.DataFrame);

# idx_pheno = findall(occursin.(pInfo.Chr, "1 2 3 4 5"));
# pInfo_subset = pInfo[idx_pheno, :];


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
K_test = Helium.readhe(joinpath(@__DIR__, "..", "test", "K_test.he"));
K_test == kinship

########
# Scan #
########
using Random
# rng = MersenneTwister(2023)
# Random.seed!(MersenneTwister(2023))
# rand(1)
# rand(Xoshiro(100))

x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] .* 1.0;
mX = reshape(x, :, 1)
mX_src = BulkLMM.transform_permute(mX; nperms = 2000, original = false);
mX_test = Helium.readhe(joinpath(@__DIR__, "..", "test", "mX_perms.he"));
mX_test == mX_src

K_eigen = BulkLMM.eigen(kinship);
eigen_test = Helium.readhe(joinpath(@__DIR__, "..", "test", "eigen_test.he"));
# rng = Xoshiro(0);shuffle(rng, x)

single_results_perms = scan(
	pheno_y,
	geno_processed,
	kinship;
	permutation_test = true,
	nperms = 2000,
);

single_results_perms2 = Helium.readhe(joinpath(@__DIR__,
	"..", "test", "scan_perms.he"));
single_results_perms2 == single_results_perms


thrs = BigRiverQTLPlots.perms_thresholds(
	single_results_perms.L_perms, [0.90, 0.95]) |> permutedims;
thrs_test = Helium.readhe(joinpath(@__DIR__,
	"..", "test", "thrs_test.he"));
thrs == thrs_test
vcat(thrs, thrs_test)

single_results = scan(
	pheno_y,
	geno_processed,
	kinship,
);


single_results3 = scan(
	pheno_y,
	geno_processed,
	kinship,
);

single_results2 = Helium.readhe(joinpath(@__DIR__,
	"..", "test", "scan_test.he"));
single_results2 == single_results.lod

########
# Plot #
########
plot_QTL(single_results.lod, gInfo)
savefig(joinpath(@__DIR__, "..", "images", "QTL_test.png"))

thr = BigRiverQTLPlots.perms_thresholds(single_results_perms.L_perms, [0.90, 0.95])
plot_QTL(single_results.lod, gInfo, thresholds = thr)
savefig(joinpath(@__DIR__, "..", "images", "QTL_thrs_test.png"))
