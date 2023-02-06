using Plots, Plots.PlotMeasures 
using BigriverPlots
using Helium

# Read data
chr_file = joinpath(@__DIR__, "..", "data", "chr.he")
pos_file = joinpath(@__DIR__, "..", "data", "pos.he")
lod_file = joinpath(@__DIR__, "..", "data", "lod.he")

vecChr = BigriverPlots.Helium.readhe(chr_file);
vecLoci = BigriverPlots.Helium.readhe(pos_file);
vecLod = BigriverPlots.Helium.readhe(lod_file);

vecSteps = BigriverPlots.get_chromosome_steps(vecLoci, vecChr)

# get unique chr id
v_chr_names = unique(vecChr)

vPos_new = BigriverPlots.get_abs_loci(vecLoci, vecChr, vecSteps)

# generate new distances coordinates

x, y = BigriverPlots.get_qtl_coord(vecLoci, vecChr, vecLod)

qtlplot(x,y, vecSteps, string.(Int.(v_chr_names)))
    
