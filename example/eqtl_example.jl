
using BulkLMM, Helium
# using CSV, DelimitedFiles, DataFrames, Statistics
# using RecipesBase
using Plots, Plots.PlotMeasures, ColorSchemes

cgrad(:blues)

########
# Data #
########
bulklmmdir = dirname(pathof(BulkLMM));

pheno_file = joinpath(bulklmmdir,"..","data/bxdData/spleen-pheno-nomissing.csv");
pheno = BulkLMM.DelimitedFiles.readdlm(pheno_file, ',', header = false);
pheno_processed = pheno[2:end, 2:(end-1)].*1.0; # exclude the header, the first (transcript ID)and the last columns (sex)

geno_file = joinpath(bulklmmdir,"..","data/bxdData/spleen-bxd-genoprob.csv")
geno = BulkLMM.DelimitedFiles.readdlm(geno_file, ',', header = false);
geno_processed = geno[2:end, 1:2:end] .* 1.0;

gmap_file = joinpath(bulklmmdir,"..","data/bxdData/gmap.csv");
gInfo = BulkLMM.CSV.read(gmap_file, BulkLMM.DataFrames.DataFrame);

phenocovar_file = joinpath(bulklmmdir,"..","data/bxdData/phenocovar.csv");
pInfo = BulkLMM.CSV.read(phenocovar_file, BulkLMM.DataFrames.DataFrame);

###########
# Kinship #
###########
kinship = calcKinship(geno_processed); 

########
# Scan #
########
multipletraits_results = scan_lite_multivar(pheno_processed, geno_processed,
                                              kinship, Threads.nthreads());
# Helium.writehe(multipletraits_results, "data/bxd/multipletraits_results.he")
# multipletraits_results = Helium.readhe("data/bxd/multipletraits_results.he");


########
# Plot #
########

using BigRiverPlots

# include(joinpath(bulklmmdir, "..", "plot_utils", "visuals_utils.jl"));

ploteQTL(multipletraits_results, pInfo, gInfo; thr = 5.0)


# pheno
# # Test
# using DataFrames
# vChr = unique(gInfo.Chr)
# combine(groupby(gInfo, :Chr), :Mb => maximum => :Steps_pheno)
# @time dfSteps = filter(:Chr => x-> [x] ⊆ vChr, pInfo)|> 
# x-> groupby(x, :Chr) |>
# x-> combine(x, :Mb => maximum => :Steps_pheno) |>
# x-> leftjoin(x, combine(groupby(gInfo, :Chr), :Mb => maximum => :Steps_geno), on =:Chr)|>
# x-> transform(x, [:Steps_pheno, :Steps_geno] => ByRow(max)=> :Steps_max)|>
# x-> transform(x, :Steps_max => cumsum => :Steps)|>
# x-> select(x, [:Chr, :Steps])

# leftjoin(gInfo, dfSteps, on= :Chr)

# sort(pInfo, :Chr, )