
using BigRiverQTLPlots
using BulkLMM, Helium
using Statistics
using Plots


##############
# BXD spleen #
##############

########
# Data #
########
bulklmmdir = dirname(pathof(BulkLMM));

pheno_file = joinpath(bulklmmdir, "..", "data", "bxdData", "spleen-pheno-nomissing.csv");
pheno = BulkLMM.DelimitedFiles.readdlm(pheno_file, ',', header = false);
pheno_processed = pheno[2:end, 2:(end-1)].*1.0; # exclude the header, the first (transcript ID)and the last columns (sex)

geno_file = joinpath(bulklmmdir,"..", "data", "bxdData", "spleen-bxd-genoprob.csv")
geno = BulkLMM.DelimitedFiles.readdlm(geno_file, ',', header = false);
geno_processed = geno[2:end, 1:2:end] .* 1.0;

gmap_file = joinpath(bulklmmdir,"..", "data", "bxdData", "gmap.csv");
gInfo = BulkLMM.CSV.read(gmap_file, BulkLMM.DataFrames.DataFrame);

phenocovar_file = joinpath(bulklmmdir,"..", "data", "bxdData", "phenocovar.csv");
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
single_results_perms = scan(
                        pheno_y, 
                        geno_processed, 
                        kinship; 
                        permutation_test = true, 
                        nperms = 2000, 
                        original = true
);                                             

single_results = scan(
                        pheno_y, 
                        geno_processed, 
                        kinship
);     
########
# Plot #
########
plot_QTL(single_results, gInfo)
savefig(joinpath(@__DIR__,"..","images","QTL_test.png"))

plot_QTL(single_results_perms, gInfo, thresholds= [0.90, 0.95])
savefig(joinpath(@__DIR__,"..","images","QTL_thrs_test.png"))


###############
# Arabidopsis #
###############

########
# Data #
########

# Read data
chr_file = joinpath(@__DIR__, "..", "data", "arabidopsisdata", "chr.he")
pos_file = joinpath(@__DIR__, "..", "data", "arabidopsisdata", "pos.he")
lod_file = joinpath(@__DIR__, "..", "data", "arabidopsisdata", "lod.he")

vecChr = BigRiverQTLPlots.Helium.readhe(chr_file);
vecLoci = BigRiverQTLPlots.Helium.readhe(pos_file);
vecLod = BigRiverQTLPlots.Helium.readhe(lod_file);

#################
# Preprocessing #
#################

vecSteps = BigRiverQTLPlots.get_chromosome_steps(vecLoci, vecChr)

# get unique chr id
v_chr_names = unique(vecChr)

vPos_new = BigRiverQTLPlots.get_pseudo_loci(vecLoci, vecChr, vecSteps)

# generate new distances coordinates

x, y = BigRiverQTLPlots.get_qtl_coord(vecLoci, vecChr, vecLod)

########
# Plot #
########

qtlplot(x,y, vecSteps, string.(Int.(v_chr_names)))
    
