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
    new_ticks = zeros(size(mytick,1))
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


function get_abs_loci(vLoc, vChr, vSteps)

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
    # sort vChr, vLoci and vLOD according to vLoci 
    mData = hcat(vChr, vLoci, vLod);
    mData = sortslices(mData, dims=1, lt=(x,y)->isless(x[2],y[2]));

    n = size(mData, 1);
    # insert Inf to obtain separation in plotting
    x = zeros(Float64, n + length(vec_chr_names));
    y = zeros(Float64, n + length(vec_chr_names));
    
    for i in eachindex(vec_chr_names)
        idx = getindex.((findall(vChr .== vec_chr_names[i])), 1);
        x[idx.+i.-1] .= vPos_new[idx]
        x[idx[end]+i] = Inf
        
        y[idx.+i.-1] .= vLod[idx]
        y[idx[end]+i] = Inf
    end
    
    x = x[1:end-1]
    y = y[1:end-1]

    return x, y
end