#=
List of the utils functions
-

- get_chromosome_steps
	Returns a vector containing the accumulated version of all maximimum loci.            

- get_pseudo_loci
	Returns a vector containing the pseudo "accumulated" loci (mB).

- get_qtl_coord
	Return coordinates vectors for plotting QTL figure.

- plot_QTL
	Generates a scatter plot for QTL analysis.        

=#

"""

function perms_thresholds(mLOD::Matrix(<: AbstractFloat), significance::Vector{<: AbstractFloat}) 

	Returns LOD thresholds.

# Arguments
- mLOD is a matrix containing the LOD scores permutations.
- significance contains the significant levels, default values are 0.10 and 0.05.

# Output
- thresholds_lod contains the LOD scores corresponding to the significant levels.

"""

function perms_thresholds(mLOD::Matrix{<: AbstractFloat}, significance::Vector = [0.10, 0.05])

	# default thresholds
	if isempty(significance)
		significance = [0.10, 0.05]
	end

	# get maximum values for each column
	max_lods = vec(mapslices(x -> maximum(x), mLOD; dims = 1))

	# get LOD values corresponding significances values
	thresholds_lod = map(x -> quantile(max_lods, x), 1 .- significance)

	return thresholds_lod
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
	vLoc_new = copy(vLoc)

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
	mData = hcat(vChr, vLoci, vLod)
	mData = sortslices(mData, dims = 1, lt = (x, y) -> isless(x[2], y[2]))

	n = size(mData, 1)
	# insert Inf to obtain separation in plotting
	x = zeros(Float64, n + length(vec_chr_names))
	y = zeros(Float64, n + length(vec_chr_names))

	vecSteps = get_chromosome_steps(vLoci, vChr)

	vLoci_new = get_pseudo_loci(vLoci, vChr, vecSteps)

	for i in eachindex(vec_chr_names)
		idx = getindex.((findall(vChr .== vec_chr_names[i])), 1)
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
get_plot_QTL_inputs(vLOD::Vector{<: AbstractFloat}, dfgInfo::DataFrame;
					chrColname::String="Chr", mbColname::String="Mb")

Obtains required inputs to generates a QTL plot.

## Arguments
- `vLOD` is the vector containing the maximum value of LOD score of each phenotype and its corresponding index.
- `dfpInfo` is a dataframe containing the phenotype information such as probeset, chromosomes names and Mb distance.
- `dfgInfo` is a dataframe containing the genotype information such as locus, cM distance, chromosomes names and Mb distance. 
- `chrColname` is the name of the column containing the chromosomes' names, default name is "Chr".
- `mbColname` is the name of the column containing the megabase DNA length, default name is "Mb". 

"""
function get_plot_QTL_inputs(vLOD::Vector{<: AbstractFloat}, dfgInfo::DataFrame;
	chrColname::String = "Chr", mbColname::String = "Mb")

	vecChr = String.(dfgInfo[:, Symbol(chrColname)])
	vecLoci = dfgInfo[:, Symbol(mbColname)]


	# get unique chr id
	v_chr_names = sortnatural(unique(vecChr))

	vecSteps = get_chromosome_steps(vecLoci, vecChr)

	x, y = get_qtl_coord(vecLoci, vecChr, vLOD)

	return x, y, vecSteps, v_chr_names
end


"""
plot_QTL(vLOD::Vector{<: AbstractFloat}, dfgInfo::DataFrame;
		chrColname::String="Chr", mbColname::String="Mb",
		thresholds::Vector{<: AbstractFloat}=[], kwargs...)

Generates a scatter plot for QTL analysis.

## Arguments
- `vLOD` is the vector containing the maximum value of LOD score.
- `dfgInfo` is a dataframe containing the genotype information such as locus, cM distance, chromosomes names and Mb distance. 
- `chrColname` is the name of the column containing the chromosomes' names, default name is "Chr".
- `mbColname` is the name of the column containing the megabase DNA length, default name is "Mb". 
- `thresholds` is <: AbstractFloat number vector containing desired LOD score thresholds for plotting.

---

plot_QTL(scanresult::NamedTuple, dfgInfo::DataFrame;
		chrColname::String="Chr", mbColname::String="Mb", 
		thresholds::Vector{<: AbstractFloat}=[], kwargs...)

Generates a scatter plot for QTL analysis.

## Arguments
- `scanresult` is NamedTuple object resulting from the `scan()` function in `BulkLMM.jl`.
- `dfgInfo` is a dataframe containing the genotype information such as locus, cM distance, chromosomes names and Mb distance. 
- `chrColname` is the name of the column containing the chromosomes' names, default name is "Chr".
- `mbColname` is the name of the column containing the megabase DNA length, default name is "Mb". 
- `significance` is <: AbstractFloat number vector containing significant levels to estimate LOD score thresholds.

If the `scanresult` does not contain a permutation matrix, the original maximum LOD scores will be plotted, and the values in 
the `significance` vector will be used as the threshold values for comparison.

"""
function plot_QTL(vLOD::Vector{<: AbstractFloat}, dfgInfo::DataFrame;
	chrColname::String = "Chr", mbColname::String = "Mb",
	thresholds::Vector = [], kwargs...)

	x, y, vecSteps, v_chr_names = get_plot_QTL_inputs(vLOD, dfgInfo;
		chrColname = chrColname, mbColname = mbColname)

	if isempty(thresholds)
		qtlplot(x, y, vecSteps, v_chr_names; kwargs...)
	else
		qtlplot(x, y, vecSteps, v_chr_names, thresholds; kwargs...)
	end

end

function plot_QTL(scanresult::NamedTuple, dfgInfo::DataFrame;
	chrColname::String = "Chr", mbColname::String = "Mb",
	significance::Vector = [], kwargs...)

	if (:L_perms in keys(scanresult))
		thrshlds = perms_thresholds(scanresult.L_perms, significance)
	else
		thrshlds = significance
	end

	plot_QTL(
		scanresult.lod, dfgInfo;
		chrColname = chrColname,
		mbColname = mbColname,
		thresholds = thrshlds,
		kwargs...,
	)

end

function plot_QTL!(vLOD::Vector{<: AbstractFloat}, dfgInfo::DataFrame;
	chrColname::String = "Chr", mbColname::String = "Mb",
	thresholds::Vector = [], kwargs...)

	x, y, vecSteps, v_chr_names = get_plot_QTL_inputs(vLOD, dfgInfo;
		chrColname = chrColname, mbColname = mbColname)

	if isempty(thresholds)
		qtlplot!(x, y, vecSteps, v_chr_names; kwargs...)
	else
		qtlplot!(x, y, vecSteps, v_chr_names, thresholds; kwargs...)
	end

end

function plot_QTL!(scanresult::NamedTuple, dfgInfo::DataFrame;
	chrColname::String = "Chr", mbColname::String = "Mb",
	thresholds::Vector = [], kwargs...)

	if (:L_perms in keys(scanresult))
		thrshlds = perms_thresholds(scanresult.L_perms, thresholds)
	else
		thrshlds = thresholds
	end

	plot_QTL!(
		scanresult.lod, dfgInfo;
		chrColname = chrColname,
		mbColname = mbColname,
		thresholds = thrshlds,
		kwargs...,
	)

end