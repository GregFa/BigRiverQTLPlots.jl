
############## 
# QTL PLOT  #
##############

"""
    Recipe for QTL plots.
"""
mutable struct QtlPlot{AbstractType}
    args::Any                                      
end

qtlplot(args...; kw...) = RecipesBase.plot(QtlPlot{typeof(args[1])}(args); kw...)

@recipe function f(h::QtlPlot) 
    # check types of the input arguments
    if length(h.args) != 4 || !(typeof(h.args[1]) <: AbstractVector) ||
        !(typeof(h.args[2]) <: AbstractVector) || !(typeof(h.args[3]) <: AbstractVector) ||
        !(typeof(h.args[4]) <: AbstractVector)  
        error("QTL Plots should be given four vectors.  Got: $(typeof(h.args))")
    end
        
    # get arguments
    x, y, steps, chr_names = h.args

    # get number of shaded area for chromosomes
    idx_bar = findall(isodd.(eachindex(steps[2:end])));

    # get maximum LOD value
    y_max = round(maximum(y[y.!=Inf]));


    # set a default value for an attribute with `-->`
    xlabel --> "Chromosome"
    ylabel --> "LOD score"
    
    marker --> 6
    markerstrokewidth --> 0.3
    
    # bottom_margin --> 0
    # right_margin --> 0
    
    guidefontsize --> 15
    fontfamily --> "Helvetica"
    
    # size --> (650, 550)
        
    # set up the subplots
    legend := false
    link := :both
    # framestyle := [:none :axes :none]
    # yaxis := false 
    xlims := (0, steps[end])
    ylims := (0, y_max)
    grid --> (:y)
    y_foreground_color_axis --> :white
    y_foreground_color_border --> :white
    x_foreground_color_border --> :white

    tickfontsize := 8
    tick_direction := :out

    xticks := (pseudoticks(steps[2:end]), chr_names)
    # yticks := (pseudotick(steps), chr_names)
    
    
  
    ############################
    # Bar separting chromosome #
    ############################
    @series begin
        seriestype := :bar
        # marker_z := lod
        # framestyle := :box
        
        label := ""
        bar_width := steps[2:end][idx_bar] .- steps[idx_bar]

        # get the seriescolor passed by the user
        alpha --> 0.2
        linecolor --> :lightsalmon#nothing
        color --> :lightsalmon
        
        
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
end