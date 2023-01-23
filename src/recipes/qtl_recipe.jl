
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

@recipe function f(h::EqtlPlot) 
    # check types of the input arguments
    if length(h.args) != 5 || !(typeof(h.args[1]) <: AbstractVector) ||
        !(typeof(h.args[2]) <: AbstractVector) || !(typeof(h.args[3]) <: AbstractVector) ||
        !(typeof(h.args[4]) <: AbstractVector) || !(typeof(h.args[5]) <: AbstractVector) 
        error("eQTL Plots should be given three vectors.  Got: $(typeof(h.args))")
    end
    # Note: if confidence or interval not symmetric, then ε should be a vector of tuple.
    # collect(zip(vCI1, vCI2));
    
    
    x, y, lod, steps, chr_names = h.args

    idx_bar = findall(isodd.(eachindex(steps[2:end])));
    y_max = round(maximum(y[y.!=Inf]));
    
    # set a default value for an attribute with `-->`
    xlabel --> "Chromosome"
    ylabel --> "LOD score"
    
    marker --> 6
    markerstrokewidth --> 0.3
    
    bottom_margin --> 0mm
    right_margin --> 0mm
    
    guidefontsize --> 15
    fontfamily --> "Helvetica"
    
    # size --> (650, 550)
        
    # set up the subplots
    legend := false
    link := :both
    # framestyle := [:none :axes :none]
    # yaxis := false 
    xlims := (0, steps[end])
    # ylims := (0, steps[end])
    grid := false
    

    tickfontsize := 8
    tick_direction := :out

    xticks := (pseudotick(steps), chr_names) 
    # yticks := (pseudotick(steps), chr_names)
    
    
    # vertical lines
    @series begin
        seriestype := :vline
        linecolor := :lightgrey
        primary := false
        # alpha := 0.5
        steps
        
    end
  
      # bar separting chromosome
      @series begin
        seriestype := :bar
        marker_z := lod
        framestyle := :box
        linecolor := :black#nothing
        
        label := ""
        bar_width:= vec_steps[2:end][idx_bar] -vec_steps[idx_bar]

        # get the seriescolor passed by the user
        alpha --> 0.2
        color -->  :lightsalmon
        
        x, y
    end  



    # # horizontal lines
    # @series begin
    #     seriestype := :hline
    #     linecolor := :lightgrey
    #     primary := false
    #     # alpha := 0.5
    #     steps
    # end
    

    # # main confidence plot
    # @series begin
    #     seriestype := :scatter
    #     marker_z := lod
    #     framestyle := :box
    #     linecolor := :black#nothing
    #     # get the seriescolor passed by the user
    #     color -->  cgrad(:blues)
    #     cbar --> true

    #     x, y
    # end  
end



#=
bar(
        pseudotick(vec_steps[2:end])[idx_bar], 
        repeat([y_max], length(idx_bar)), 
        color = :lightsalmon, 
        lc = :lightsalmon, 
        # fill_z = 4:-1:1, 
        alpha = 0.2, 
        label = "", 
        bar_width = vec_steps[2:end][idx_bar] -vec_steps[idx_bar],
    )

plot!(
    x, y, 
    grid = (:y),
    legend= false,
    xticks = (pseudotick(vec_steps[2:end]), string.(Int.(vec_chr_names))),
    c = :skyblue4,
    xlim = (0, maximum(x[x.!=Inf])),
    ylim= (0, y_max), 
    )
vline!(
    vec_steps[2:end],
    linecolor= :lightgrey
)
=#