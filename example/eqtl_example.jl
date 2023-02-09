
using BulkLMM, Helium
# using CSV, DelimitedFiles, DataFrames, Statistics
using RecipesBase, Plots, Plots.PlotMeasures, ColorSchemes

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
multipletraits_results = scan_lite_multivar(pheno_processed, geno_processed, kinship, Threads.nthreads());
Helium.writehe(multipletraits_results, "data/bxd/multipletraits_results.he")


include(joinpath(bulklmmdir, "..", "plot_utils", "visuals_utils.jl"));
plot_eQTL(multipletraits_results, pheno, gInfo, pInfo; thr = 5.0)



bulklmmdir = dirname(pathof(BulkLMM));
pheno_file = joinpath(bulklmmdir,"..","data/bxdData/spleen-pheno-nomissing.csv");
pheno = readdlm(pheno_file, ',', header = false);
pheno_processed = pheno[2:end, 2:(end-1)].*1.0; # exclude the header, the first (transcript ID)and the last columns (sex)



@time single_results = scan(pheno_y, geno_processed, kinship);



# using Plots, Plots.PlotMeasures, ColorSchemes 
# using BigRiverPlots
# using Helium

# # Read data
# chr_file = joinpath(@__DIR__, "..", "data", "arabidopsisdata", "chr.he")
# pos_file = joinpath(@__DIR__, "..", "data", "arabidopsisdata", "pos.he")
# lod_file = joinpath(@__DIR__, "..", "data", "arabidopsisdata", "lod.he")

# vecChr = BigRiverPlots.Helium.readhe(chr_file);
# vecLoci = BigRiverPlots.Helium.readhe(pos_file);
# vecLod = BigRiverPlots.Helium.readhe(lod_file);

# vecSteps = BigRiverPlots.get_chromosome_steps(vecLoci, vecChr)

# # get unique chr id
# v_chr_names = unique(vecChr)

# vPos_new = BigRiverPlots.get_abs_loci(vecLoci, vecChr, vecSteps)

# # generate new distances coordinates

# x, y = BigRiverPlots.get_qtl_coord(vecLoci, vecChr, vecLod)

# qtlplot(x,y, vecSteps, string.(Int.(v_chr_names)))
    



# using RecipesBase, Plots, Plots.PlotMeasures, ColorSchemes
# include(joinpath(bulklmmdir, "..", "plot_utils", "visuals_utils.jl"));

# multiple_results_allTraits = scan_lite_multivar(pheno_processed, geno_processed, kinship, Threads.nthreads());
# gmap_file = joinpath(bulklmmdir,"..","data/bxdData/gmap.csv");
# gInfo = CSV.read(gmap_file, DataFrame);
# phenocovar_file = joinpath(bulklmmdir,"..","data/bxdData/phenocovar.csv");
# pInfo = CSV.read(phenocovar_file, DataFrame);

# plot_eQTL(multiple_results_allTraits, pheno, gInfo, pInfo; thr = 5.0)
