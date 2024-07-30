"""
plot_manhattan(vLOD::Vector{<: AbstractFloat}, dfgInfo::DataFrame;
		chrColname::String="Chr", posColname::String="Mb",
		thresholds::Vector{<: AbstractFloat}=[], kwargs...)

Generates a scatter plot for QTL analysis.

## Arguments
- `vLOD` is the vector containing the maximum value of LOD score.
- `dfgInfo` is a dataframe containing the genotype information such as locus, cM distance, chromosomes names and Mb distance. 
- `chrColname` is the name of the column containing the chromosomes' names, default name is "Chr".
- `posColname` is the name of the column containing the megabase DNA length, default name is "Mb". 
- `thresholds` is <: AbstractFloat number vector containing desired LOD score thresholds for plotting.

---

plot_manhattan(scanresult::NamedTuple, dfgInfo::DataFrame;
		chrColname::String="Chr", posColname::String="Mb", 
		thresholds::Vector{<: AbstractFloat}=[], kwargs...)

Generates a scatter plot for QTL analysis.

## Arguments
- `scanresult` is NamedTuple object resulting from the `scan()` function in `BulkLMM.jl`.
- `dfgInfo` is a dataframe containing the genotype information such as locus, cM distance, chromosomes names and Mb distance. 
- `chrColname` is the name of the column containing the chromosomes' names, default name is "Chr".
- `posColname` is the name of the column containing the megabase DNA length, default name is "Mb". 
- `significance` is <: AbstractFloat number vector containing significant levels to estimate LOD score thresholds.

If the `scanresult` does not contain a permutation matrix, the original maximum LOD scores will be plotted, and the values in 
the `significance` vector will be used as the threshold values for comparison.

"""
function plot_manhattan(vLOD::Vector{<: AbstractFloat}, dfgInfo::DataFrame;
	chrColname::String = "Chr", posColname::String = "Mb",
	thresholds::Vector = [], kwargs...)

	x, y, vecSteps, v_chr_names = get_plot_QTL_inputs(vLOD, dfgInfo;
		chrColname = chrColname, posColname = posColname)

	if isempty(thresholds)
		manhattanplot(x, y, vecSteps, v_chr_names; kwargs...)
	else
		manhattanplot(x, y, vecSteps, v_chr_names, thresholds; kwargs...)
	end

end

function plot_manhattan(scanresult::NamedTuple, dfgInfo::DataFrame;
	chrColname::String = "Chr", posColname::String = "Mb",
	significance::Vector = [], kwargs...)

	if (:L_perms in keys(scanresult))
		thrshlds = perms_thresholds(scanresult.L_perms, significance)
	else
		thrshlds = significance
	end

	plot_manhattan(
		scanresult.lod, dfgInfo;
		chrColname = chrColname,
		posColname = posColname,
		thresholds = thrshlds,
		kwargs...,
	)

end

function plot_manhattan!(vLOD::Vector{<: AbstractFloat}, dfgInfo::DataFrame;
	chrColname::String = "Chr", posColname::String = "Mb",
	thresholds::Vector = [], kwargs...)

	x, y, vecSteps, v_chr_names = get_plot_QTL_inputs(vLOD, dfgInfo;
		chrColname = chrColname, posColname = posColname)

	if isempty(thresholds)
		manhattanplot!(x, y, vecSteps, v_chr_names; kwargs...)
	else
		manhattanplot!(x, y, vecSteps, v_chr_names, thresholds; kwargs...)
	end

end

function plot_manhattan!(scanresult::NamedTuple, dfgInfo::DataFrame;
	chrColname::String = "Chr", posColname::String = "Mb",
	thresholds::Vector = [], kwargs...)

	if (:L_perms in keys(scanresult))
		thrshlds = perms_thresholds(scanresult.L_perms, thresholds)
	else
		thrshlds = thresholds
	end

	plot_manhattan!(
		scanresult.lod, dfgInfo;
		chrColname = chrColname,
		posColname = posColname,
		thresholds = thrshlds,
		kwargs...,
	)

end

