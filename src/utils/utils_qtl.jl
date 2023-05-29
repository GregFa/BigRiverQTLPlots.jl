#=
List of the utils functions
- get_chromosome_steps
    Returns a vector containing the accumulated version of all maximimum loci.            

- get_pseudo_loci
    Returns a vector containing the pseudo "accumulated" loci (mB).

- get_qtl_coord
    Return coordinates vectors for plotting QTL figure.
=#


"""

get_chromosome_steps(vLoc, vChr) => Vector(::Float64)

Returns a vector containing the accumulated version of 
all maximimum loci.

# Arguments
- vLoc contains the loci 
- vChr contains the chromosome names

"""
function get_chromosome_steps(vLoc, vChr)
    
    # get unique chromosome name
    vec_chr_names = unique(vChr)

    # initiate steps vector
    vec_steps = zeros(length(vec_chr_names))

    for i in eachindex(vec_chr_names) 
        vec_steps[i] = maximum(vLoc[findall(vChr .== vec_chr_names[i])]) 
    end
    
    vec_steps = vcat([0], accumulate(+, vec_steps))

    return vec_steps
end
"""
get_pseudo_loci(vLoc, vChr, vSteps) => Vector(::Float64)

Returns a vector containing the pseudo "accumulated" loci (mB).

# Arguments
- vLoc contains the loci 
- vChr contains the chromosome names
- vSteps contains the accumulated version of all maximimum loci.

"""
function get_pseudo_loci(vLoc, vChr, vSteps)

    # get unique chromosome name
    vec_chr_names = unique(vChr)

    # generate new distances coordinates
    vLoc_new = copy(vLoc);
    
    for i in eachindex(vec_chr_names)
        idx = findall(vChr .== vec_chr_names[i])
        vLoc_new[idx] = vLoc_new[idx] .+ vSteps[i]
    end

    return vLoc_new
end

"""
get_qtl_coord(vLoci, vChr, vLod)

Return coordinates vectors for plotting QTL figure.

# Arguments
- vLoc contains the loci 
- vChr contains the chromosome names
- vLod contains the LOD scores


"""
function get_qtl_coord(vLoci, vChr, vLod)
    # get unique chromosome name
    vec_chr_names = unique(vChr)

    # sort vChr, vLoci and vLOD according to vLoci 
    mData = hcat(vChr, vLoci, vLod);
    mData = sortslices(mData, dims=1, lt=(x,y)->isless(x[2],y[2]));

    n = size(mData, 1);
    # insert Inf to obtain separation in plotting
    x = zeros(Float64, n + length(vec_chr_names));
    y = zeros(Float64, n + length(vec_chr_names));
    
    vecSteps = get_chromosome_steps(vLoci, vChr)

    vLoci_new = get_pseudo_loci(vLoci, vChr, vecSteps)

    for i in eachindex(vec_chr_names)
        idx = getindex.((findall(vChr .== vec_chr_names[i])), 1);
        x[idx.+i.-1] .= vLoci_new[idx]
        x[idx[end]+i] = Inf
        
        y[idx.+i.-1] .= vLod[idx]
        y[idx[end]+i] = Inf
    end
    
    x = x[1:end-1]
    y = y[1:end-1]

    return x, y
end

"""
get_plotqtl_inputs(mLOD::Array{Float64, 2}, dfgInfo::DataFrame;
                    chrColname::String = "Chr", mbColname::String = "Mb")

Generates a scatter plot for QTL analysis.

## Arguments
- `mLOD` is the matrix containing the maximum value of LOD score of each phenotype and its corresponding index.
- `dfpInfo` is a dataframe containing the phenotype information such as probeset, chromosomes names and Mb distance.
- `dfgInfo` is a dataframe containing the genotype information such as locus, cM distance, chromosomes names and Mb distance. 
- `chrColname` is the name of the column containing the chromosomes' names, default name is "Chr".
- `mbColname` is the name of the column containing the megabase DNA length, default name is "Mb". 

"""
function get_plotqtl_inputs(mLOD::Array{Float64, 2}, dfgInfo::DataFrame;
                            chrColname::String = "Chr", mbColname::String = "Mb")

    vecChr = String.(dfgInfo[:, Symbol(chrColname)]);
    vecLoci = dfgInfo[:, Symbol(mbColname)];
    vecLod = mLOD[:,1];

    # get unique chr id
    v_chr_names = sortnatural(unique(vecChr));

    vecSteps = get_chromosome_steps(vecLoci, vecChr);

    x, y = get_qtl_coord(vecLoci, vecChr, vecLod);

    return x, y, vecSteps, v_chr_names
end


"""
plotQTL(mLOD::Array{Float64, 2}, dfgInfo::DataFrame;
        chrColname::String = "Chr", mbColname::String = "Mb", 
        thresholds = [],
        kwargs...)

Generates a scatter plot for QTL analysis.

## Arguments
- `mLOD` is the matrix containing the maximum value of LOD score of each phenotype and its corresponding index.
- `dfpInfo` is a dataframe containing the phenotype information such as probeset, chromosomes names and Mb distance.
- `dfgInfo` is a dataframe containing the genotype information such as locus, cM distance, chromosomes names and Mb distance. 
- `chrColname` is the name of the column containing the chromosomes' names, default name is "Chr".
- `mbColname` is the name of the column containing the megabase DNA length, default name is "Mb". 
- `thresholds` is real number vector containing desired quantile thresholds for plotting.

"""
function plotQTL(mLOD::Array{Float64, 2}, dfgInfo::DataFrame;
                chrColname::String = "Chr", mbColname::String = "Mb", 
                thresholds = [],
                kwargs...)

    x, y, vecSteps, v_chr_names = get_plotqtl_inputs(mLOD, dfgInfo;
                                    chrColname = chrColname, mbColname = mbColname)

                
    if !isempty(thresholds)
        max_lods = vec(mapslices(x -> maximum(x), mLOD; dims = 1));
        thrs = map(x -> quantile(max_lods, x), thresholds);
        qtlplot(x,y, vecSteps, v_chr_names, thrs, kwargs...)
    else
        qtlplot(x,y, vecSteps, v_chr_names, kwargs...)
    end

end
