using Helium
using BigRiverQTLPlots

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
# Save #
########

BigRiverQTLPlots.Helium.writehe(
    reshape(x, length(x), 1), 
    "data/loci.he";
)
BigRiverQTLPlots.Helium.writehe(
    reshape(y, length(y), 1),
     "data/lod.he"
) 

    
