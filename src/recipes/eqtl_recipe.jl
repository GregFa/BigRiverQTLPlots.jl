
############## 
# EQTL PLOT  #
##############

"""
    Recipe for eQTL plots.
"""
# mutable struct EqtlPlot{AbstractType}
#     args::Any                                      
# end

# eqtlplot(args...; kw...) = RecipesBase.plot(EqtlPlot{typeof(args[1])}(args); kw...)

@userplot EQTLPlot

@recipe function f(h::EQTLPlot) 
    # check types of the input arguments
    if length(h.args) != 5 || !(typeof(h.args[1]) <: AbstractVector) ||
        !(typeof(h.args[2]) <: AbstractVector) || !(typeof(h.args[3]) <: AbstractVector) ||
        !(typeof(h.args[4]) <: AbstractVector) || !(typeof(h.args[5]) <: AbstractVector) 
        error("eQTL Plots should be given three vectors.  Got: $(typeof(h.args))")
    end
    # Note: if confidence or interval not symmetric, then ε should be a vector of tuple.
    # collect(zip(vCI1, vCI2));
    
    
    x, y, lod, steps, chr_names = h.args
    
    # set a default value for an attribute with `-->`
    xlabel --> "eQTL Position (Chromosome)"
    ylabel --> "Transcript Position (Chromosome)"
    
    marker --> 6
    markerstrokewidth --> 0.3
    
    bottom_margin --> (0, :mm)
    right_margin --> (0, :mm)
    
    guidefontsize --> 15
    fontfamily --> "Helvetica"
    
    size --> (650, 550)
        
    # set up the subplots
    legend := false
    link := :both
    # framestyle := [:none :axes :none]
    # yaxis := false 
    xlims := (0, steps[end])
    ylims := (0, steps[end])
    grid := false
    

    tickfontsize := 8
    tick_direction := :out

    xticks := (pseudoticks(steps), chr_names) 
    yticks := (pseudoticks(steps), chr_names)
    
    
    # vertical lines
    @series begin
        seriestype := :vline
        linecolor := :lightgrey
        primary := false
        # alpha := 0.5
        steps
        
    end
    
    # horizontal lines
    @series begin
        seriestype := :hline
        linecolor := :lightgrey
        primary := false
        # alpha := 0.5
        steps
    end
    

    # main scatter plot
    @series begin
        seriestype := :scatter
        marker_z := lod
        framestyle := :box
        linecolor := :black#nothing
        # get the seriescolor passed by the user
        color -->  :blues #cgrad(:blues)
        cbar --> true

        x, y
    end  
end
