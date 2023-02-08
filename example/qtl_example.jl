using Plots, Plots.PlotMeasures 
using BigRiverPlots
using Helium

# Read data
chr_file = joinpath(@__DIR__, "..", "data", "arabidopsisdata", "chr.he")
pos_file = joinpath(@__DIR__, "..", "data", "arabidopsisdata", "pos.he")
lod_file = joinpath(@__DIR__, "..", "data", "arabidopsisdata", "lod.he")

vecChr = BigRiverPlots.Helium.readhe(chr_file);
vecLoci = BigRiverPlots.Helium.readhe(pos_file);
vecLod = BigRiverPlots.Helium.readhe(lod_file);

vecSteps = BigRiverPlots.get_chromosome_steps(vecLoci, vecChr)

# get unique chr id
v_chr_names = unique(vecChr)

vPos_new = BigRiverPlots.get_abs_loci(vecLoci, vecChr, vecSteps)

# generate new distances coordinates

x, y = BigRiverPlots.get_qtl_coord(vecLoci, vecChr, vecLod)

qtlplot(x,y, vecSteps, string.(Int.(v_chr_names)))
    
