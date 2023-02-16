#=
List of the utils functions
- pseudoticks
    Returns coordinates of the new ticks.

- sortnatural
    Natural sort a string vector accounting for numeric sorting.

=#


"""
**pseudoticks** -*Function*.
    pseudotick(mytick::Vector{Float64}) => Vector{Float64}
Returns coordinates of the new ticks. 
"""
function pseudoticks(myticks::Vector{Float64})
    new_ticks = zeros(size(myticks,1))
    for i in 1:size(myticks,1)
        if i == 1
            new_ticks[i] = myticks[i]/2
        else
            new_ticks[i] = myticks[i-1] + (myticks[i] - myticks[i-1])/2
        end 
    end   
    return new_ticks
end

"""

sortnatural(x::Vector{String}) => Vector(::String)

Natural sort a string vector accounting for numeric sorting.

## Example
```julia
julia> myvec = ["12", "2", "1", "Y", "10", "X"];
julia> sort(myvec)
6-element Vector{String}:
 "1"
 "10"
 "12"
 "2"
 "X"
 "Y"
 julia> sortnatural(myvec)
 6-element Vector{String}:
  "1"
  "2"
  "10"
  "12"
  "X"
  "Y"
```
"""
function sortnatural(x::Vector{String})
    f = text -> all(isnumeric, text) ? Char(parse(Int, text)) : text
    sorter = key -> join(f(m.match) for m in eachmatch(r"[0-9]+|[^0-9]+", key))
    sort(x, by=sorter)
end




# """

# get_chromosome_steps(vLoc, vChr) => Vector(::Float64)

# Returns a vector containing the accumulated version of 
# all maximimum loci.

# # Arguments
# - vLoc contains the loci 
# - vChr contains the chromosome names

# """
# function get_chromosome_steps(vLoc, vChr)
    
#     # get unique chromosome name
#     vec_chr_names = unique(vChr)

#     # initiate steps vector
#     vec_steps = zeros(length(vec_chr_names))

#     for i in eachindex(vec_chr_names) 
#         vec_steps[i] = maximum(vLoc[findall(vChr .== vec_chr_names[i])]) 
#     end
    
#     vec_steps = vcat([0], accumulate(+, vec_steps))

#     return vec_steps
# end
# """
# get_pseudo_loci(vLoc, vChr, vSteps) => Vector(::Float64)

# Returns a vector containing the pseudo "accumulated" loci (mB).

# # Arguments
# - vLoc contains the loci 
# - vChr contains the chromosome names
# - vSteps contains the accumulated version of all maximimum loci.

# """
# function get_pseudo_loci(vLoc, vChr, vSteps)

#     # get unique chromosome name
#     vec_chr_names = unique(vChr)

#     # generate new distances coordinates
#     vLoc_new = copy(vLoc);
    
#     for i in eachindex(vec_chr_names)
#         idx = findall(vChr .== vec_chr_names[i])
#         vLoc_new[idx] = vLoc_new[idx] .+ vSteps[i]
#     end

#     return vLoc_new
# end

# """
# get_qtl_coord(vLoci, vChr, vLod)

# Return coordinates vectors for plotting QTL figure.

# # Arguments
# - vLoc contains the loci 
# - vChr contains the chromosome names
# - vLod contains the LOD scores


# """
# function get_qtl_coord(vLoci, vChr, vLod)
#     # get unique chromosome name
#     vec_chr_names = unique(vChr)

#     # sort vChr, vLoci and vLOD according to vLoci 
#     mData = hcat(vChr, vLoci, vLod);
#     mData = sortslices(mData, dims=1, lt=(x,y)->isless(x[2],y[2]));

#     n = size(mData, 1);
#     # insert Inf to obtain separation in plotting
#     x = zeros(Float64, n + length(vec_chr_names));
#     y = zeros(Float64, n + length(vec_chr_names));
    
#     vecSteps = get_chromosome_steps(vLoci, vChr)

#     vLoci_new = get_pseudo_loci(vLoci, vChr, vecSteps)

#     for i in eachindex(vec_chr_names)
#         idx = getindex.((findall(vChr .== vec_chr_names[i])), 1);
#         x[idx.+i.-1] .= vLoci_new[idx]
#         x[idx[end]+i] = Inf
        
#         y[idx.+i.-1] .= vLod[idx]
#         y[idx[end]+i] = Inf
#     end
    
#     x = x[1:end-1]
#     y = y[1:end-1]

#     return x, y
# end


# """
# get_eQTL_accMb(mlodmax::Matrix{Float64}, dfpInfo::DataFrame, dfgInfo::DataFrame; 
#                chrColname::String = "Chr", mbColname::String = "Mb", thrs::Float64 = 5.0)=> 


# Returns pseudo "acculumated" loci for plotting, according to the phenotype and genotype annotations.

# ## Arguments
# - `mlodmax` is the matrix containing the maximum value of LOD score of each phenotype and its corresponding index
# - `dfpInfo` is a dataframe containing the phenotype information such as probeset, chromosomes names and Mb distance
# - `dfgInfo` is a dataframe containing the genotype information such as locus, cM distance, chromosomes names and Mb distance  
# - `chrColname` column name containing the chromosomes information, default is `"Chr"`, in the dataframes. 
# - `mbColname` column name containing the Mb distance information, default is `"Mb"`, in the dataframes.
# - `thrs` is the LOD threshold value, default is `5.0``.

# ## Output
# - `pheno_gmap_lod.acc_geno_mb` vector contains pseudo "accumulated" loci (mB) based on genotype.
# - `pheno_gmap_lod.acc_phenocovar_mb` vector contains pseudo "accumulated" loci (mB) based on phenotype.
# - `pheno_gmap_lod.maxlod` vector contains the maximum LOD scores above the threshold value `thrs`.  
# - `steps` contains the accumulated version of all maximimum loci above the threshold value `thrs`.
# - `vChrNames` contains the chromosome names.


# """
# function get_eQTL_accMb(mlodmax::Matrix{Float64}, dfpInfo::DataFrame, dfgInfo::DataFrame; 
#                         chrColname::String = "Chr", mbColname::String = "Mb", thr::Float64 = 5.0)
#     # match chromosomes in pheno dataframe according to chromosomes list in geno dataframe
#     dfpInfo_filtered = copy(dfpInfo) #match_chrs_pheno_to_geno(dfpInfo, dfgInfo)

#     # prepare filtered pheno dataframe to compute accumulated mb distance for plotting
#     rename!(dfpInfo_filtered, Dict(Symbol(mbColname) => :phenocovar_mb));
#     rename!(dfpInfo_filtered, Dict(Symbol(chrColname) => :phenocovar_chr)); 
#     dfpInfo_filtered[:, :acc_phenocovar_mb] .= dfpInfo_filtered.phenocovar_mb;

#     # get a copy of genotype info dataframe
#     gmap = copy(dfgInfo) 
#     gmap.acc_geno_mb .= gmap[!, mbColname]
#     rename!(gmap, Dict(Symbol(mbColname) => :geno_mb));
#     rename!(gmap, Dict(Symbol(chrColname) => :geno_chr)); 
    
#     # get unique chromosomes names from genotype info dataframe
#     vChrNames = unique(gmap.geno_chr);
#     vChrNames = sortnatural(String.(vChrNames))

#     # initiate steps vector
#     steps = Array{Float64}(undef, length(vChrNames))

#     # if more than one chromosome
#     if length(vChrNames) > 1 
#         # compute accumulated mb distance
#         for i in 1:length(vChrNames)-1 
            
#             # get temp matrix based on a chromosome name
#             phenotemp = filter(:phenocovar_chr => x-> x == vChrNames[i], dfpInfo_filtered)
#             genotemp = filter(:geno_chr => x-> x == vChrNames[i], gmap)

#             # get maximum distance for this chromosome
#             steps[i] = max(
#                         maximum(genotemp.acc_geno_mb),
#                         maximum(phenotemp.acc_phenocovar_mb) 
#                     )

#             # calculate the accumulated distance for the next chromosome
#             nextchr_phenotemp = view(dfpInfo_filtered, dfpInfo_filtered.phenocovar_chr .== vChrNames[i+1], :)
#             nextchr_genotemp = view(gmap, gmap.geno_chr .== vChrNames[i+1], :) 

#             nextchr_phenotemp.acc_phenocovar_mb .= nextchr_phenotemp.phenocovar_mb .+ steps[i]
#             nextchr_genotemp.acc_geno_mb .= nextchr_genotemp.geno_mb .+ steps[i]
#         end
#     end

#     steps[end] = max(
#                   maximum(dfpInfo_filtered.acc_phenocovar_mb), 
#                   maximum(gmap.acc_geno_mb)
#                 )

#     # get the corresponding index of the maximum values of LOD score of each phenotype
#     idxlodmax = trunc.(Int, mlodmax[:,1])
#     # concatenate results and gmap info
#     pheno_gmap_lod = hcat(dfpInfo_filtered, gmap[idxlodmax,:], DataFrame(idx = idxlodmax, maxlod = mlodmax[:,2]));

#     # filter results according to the chromosome names in geno file
#     pheno_gmap_lod = filter(:phenocovar_chr => in(vChrNames), pheno_gmap_lod);

#     # filter according to LOD threshold
#     pheno_gmap_lod = filter(row -> row.maxlod > thr, pheno_gmap_lod);

#     return  pheno_gmap_lod.acc_geno_mb, pheno_gmap_lod.acc_phenocovar_mb,  pheno_gmap_lod.maxlod, steps, vChrNames
# end



# function eQTLplot(multiLODs::Array{Float64, 2}, gmap::DataFrame, phenocovar::DataFrame;
#     thr::Float64 = 5.0, kwargs...)

# maxLODs_allTraits = mapslices(x -> findmax(x), multiLODs; dims = 1);
# maxLODs_allTraits = reduce(vcat, vec(maxLODs_allTraits));
# lodc = Array{Float64, 2}(undef, size(multiLODs, 2), 2);

# for i in 1:size(multiLODs, 2)
# lodc[i, 1] = maxLODs_allTraits[i][2];
# lodc[i, 2] = maxLODs_allTraits[i][1];
# end

# x, y, z, mysteps, mychr = get_eQTL_accMb(
#                     lodc, 
#                     phenocovar,
#                     gmap;
#                     thr = thr,
#                     kwargs...
#                   )

# eqtlplot(x, y, z, mysteps, mychr, kwargs...)

# end
