
############## 
# QTL PLOT  #
##############

"""
	Recipe for QTL plots.
"""
# mutable struct QtlPlot{AbstractType}
#     args::Any                                      
# end

# qtlplot(args...; kw...) = RecipesBase.plot(QtlPlot{typeof(args[1])}(args); kw...)

@userplot QTLPlot

@recipe function f(h::QTLPlot;
	barcolor = :lightsalmon)
	# check types of the input arguments
	if length(h.args) < 4 || !(typeof(h.args[1]) <: AbstractVector) ||
	   !(typeof(h.args[2]) <: AbstractVector) || !(typeof(h.args[3]) <: AbstractVector) ||
	   !(typeof(h.args[4]) <: AbstractVector)
		error("QTL Plots should be given at least four vectors.  Got: $(typeof(h.args))")
	end

    #############
	# Arguments #
	#############
	# get arguments
	if length(h.args) == 4
		x, y, steps, chr_names = h.args
	else
		x, y, steps, chr_names, thresh = h.args
	end

	#################
	# Bars location #
	#################
	
    # get number of shaded area for chromosomes
	idx_bar = findall(isodd.(eachindex(steps[2:end])))

	###################
	# Axis attributes #
	###################

	# get maximum LOD value
	if length(h.args) == 4
		y_max = 1.25 * round(maximum(y[y.!=Inf]))
	else
		y_max = 1.25 * round(maximum(vcat(y[y.!=Inf], thresh)))
	end

	# set a default value for an attribute with `-->`
	xlabel --> "Locus (Chromosome)"
	ylabel --> "LOD score"

	marker --> 6
	markerstrokewidth --> 0.3

	# bottom_margin --> (0, :mm)
	right_margin --> (3, :mm)

	guidefontsize --> 15
	fontfamily --> "Helvetica"

	# size --> (650, 550)

	# set up the subplots
	legend --> false
	link := :both
	# framestyle := [:none :axes :none]
	# yaxis := false 
	xlims --> (0, steps[end])
	ylims --> (0, y_max)
	grid --> (:y)
	y_foreground_color_axis --> :white
	y_foreground_color_border --> :white
	x_foreground_color_border --> :white

	tickfontsize := 8
	tick_direction := :out

	xticks := (pseudoticks(steps[2:end]), chr_names)
	# yticks := (pseudotick(steps), chr_names)


	#############################
	# Bar separating chromosome #
	#############################
	@series begin
		seriestype := :bar
		# marker_z := lod
		# framestyle := :box

		label := ""
		bar_width := steps[2:end][idx_bar] .- steps[idx_bar]

		# get the seriescolor passed by the user
		alpha --> 0.2
		linecolor --> barcolor#nothing
		color --> barcolor


		pseudoticks(steps[2:end])[idx_bar], repeat([y_max], length(idx_bar))
	end

	##############
	# LOD values #
	##############
	@series begin
		seriestype := :path
		markershape := :none
		# color --> :skyblue4
		linecolor --> :skyblue4

		x, y

	end

	##################
	# Vertical lines #
	##################
	@series begin
		seriestype := :vline
		linecolor := :lightgrey
		primary := false
		# alpha := 0.5
		steps[1:end]

	end

	####################
	# Horizontal lines #
	####################
	if length(h.args) == 5
		@series begin
			seriestype := :hline
			linecolor --> :red
			linestyle --> :dash
			primary := false
			# alpha := 0.5
			thresh
		end
	end
end

################### 
# MANHATTAN PLOT  #
###################

"""
	Recipe for QTL plots.
"""
# mutable struct QtlPlot{AbstractType}
#     args::Any                                      
# end

# qtlplot(args...; kw...) = RecipesBase.plot(QtlPlot{typeof(args[1])}(args); kw...)

@userplot ManhattanPlot

@recipe function f(h::ManhattanPlot;
	manhattancolor = ["#bcbddc", "#756bb1"])
	# check types of the input arguments
	if length(h.args) < 4 || !(typeof(h.args[1]) <: AbstractVector) ||
	   !(typeof(h.args[2]) <: AbstractVector) || !(typeof(h.args[3]) <: AbstractVector) ||
	   !(typeof(h.args[4]) <: AbstractVector)
		error("QTL Plots should be given at least four vectors.  Got: $(typeof(h.args))")
	end

	#############
	# Arguments #
	#############
	if length(h.args) == 4
		x, y, steps, chr_names = h.args
	else
		x, y, steps, chr_names, thresh = h.args
	end


	#######################
	# Binary color style  #
	#######################
	# we use 2 alternating colors, for all chromosomes

	# ininitialize vector indicating chromosome color
	n = length(y)
	is_chr_odd = repeat([true], n)

	# find indices chromosme frontier
	idx_Inf = vcat([1], findall(y .== Inf))

	# assign true for odd chromosome position and false otherwise
	for i in 1:length(idx_Inf)-1
		is_chr_odd[idx_Inf[i]:idx_Inf[i+1]] .= isodd(i)
	end

	###################
	# Axis attributes #
	###################

	# get maximum LOD value
	if length(h.args) == 4
		y_max = 1.25 * round(maximum(y[y.!=Inf]))
	else
		y_max = 1.25 * round(maximum(vcat(y[y.!=Inf], thresh)))
	end

	# set a default value for an attribute with `-->`
	xlabel --> "Locus (Chromosome)"
	ylabel --> "LOD score"


	# bottom_margin --> (0, :mm)
	right_margin --> (3, :mm)

	guidefontsize --> 15
	fontfamily --> "Helvetica"

	# size --> (650, 550)

	# set up the subplots
	legend --> false
	# link := :both
	# framestyle := [:none :axes :none]
	# yaxis := false 
	xlims --> (0, steps[end])
	ylims --> (0, y_max)
	grid --> (:y)
	y_foreground_color_axis --> :white
	y_foreground_color_border --> :white
	x_foreground_color_border --> :white

	tickfontsize := 8
	tick_direction := :out

	xticks := (pseudoticks(steps[2:end]), chr_names)
	# yticks := (pseudotick(steps), chr_names)


	##############
	# LOD values #
	##############
	@series begin
		seriestype := :scatter
		markershape --> :circle
		markercolor --> ifelse.(is_chr_odd, manhattancolor[1], manhattancolor[2])
		markersize --> 3
		markerstrokewidth --> 0.1

		x, y

	end


	####################
	# Horizontal lines #
	####################
	if length(h.args) == 5
		@series begin
			seriestype := :hline
			linecolor --> :red
			linestyle --> :dash
			primary := false
			# alpha := 0.5
			thresh
		end
	end
end
