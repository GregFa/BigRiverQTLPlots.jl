
using BigRiverQTLPlots
using FlxQTL
using Random, Statistics
using Plots
using DelimitedFiles, DataFrames, Helium

###############
# Arabidopsis #
###############

########
# Data #
########
flxqtldir = dirname(pathof(FlxQTL));

pheno_file = joinpath(flxqtldir, "..", "data", "Arabidopsis_fitness.csv");
pheno = readdlm(pheno_file,',';skipstart=1); # skip to read the first row (column names) to obtain a matrix only

geno_file = joinpath(flxqtldir, "..", "data", "Arabidopsis_genotypes.csv");
geno = readdlm(geno_file,',';skipstart=1); 

markerinfo_file = joinpath(flxqtldir, "..", "data", "Arabidopsis_markerinfo_1d.csv");
markerinfo = readdlm(markerinfo_file,',';skipstart=1);

gInfo = gInfo = DataFrame(
    Locus = markerinfo[:, 1];
    Chr = string.(markerinfo[:, 2]), 
    cM = markerinfo[:,3], 
);

#################
# Preprocessing #
#################
Y=convert(Array{Float64,2},pheno'); #convert from transposed one to a Float64 matrix
Ystd=(Y.-mean(Y,dims=2))./std(Y,dims=2); # sitewise normalization

XX=Markers(markerinfo[:,1],markerinfo[:,2],markerinfo[:,3],geno') # marker names, chromosomes, marker positions, genotypes

Z=hcat(ones(6),vcat(-ones(3),ones(3)))
m,q = size(Z)

###########
# Kinship #
###########
K = shrinkg(kinshipMan,10,XX.X);
T,λ = K2eig(K);

########
# Scan #
########

LODs,B,est0 = geneScan(1,T,λ,Ystd,XX); # MLMM

thrs = BigRiverQTLPlots.perms_thresholds(single_results_perms.L_perms, [0.10, 0.05]);
Helium.writehe(reshape(thrs, :, 1), joinpath(@__DIR__, "..", "test", "data", "thresholds.he"))
########
# Plot #
########

plot_QTL(Float64.(LODs), gInfo; mbColname = "cM")


plot_QTL(single_results.lod, gInfo)
# savefig(joinpath(@__DIR__, "..", "images", "QTL_example.png"))

plot_QTL(single_results_perms, gInfo, significance = [0.10, 0.05])
savefig(joinpath(@__DIR__, "..", "images", "QTL_thrs_example.png"))
savefig(joinpath(@__DIR__, "..", "images", "QTL_thrs_example.svg"))

