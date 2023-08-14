module BigRiverQTLPlots

    # dependent packages 
    using DataFrames, Statistics
    using RecipesBase
    
    # utils functions
    include("utils.jl");
    export sortnatural, pseudoticks,calcKinship2
    

    # qtl functions
    include("./qtl/qtl_recipe.jl");
    export qtlplot, QTLPlot, manhattanplot, ManhattanPlot    

    include("./qtl/qtl_helpers.jl");
    export  perms_thresholds, get_plot_QTL_inputs
    
    include("./qtl/plot_QTL.jl");
    export plot_QTL, plot_QTL!

    include("./qtl/plot_manhattan.jl");
    export plot_manhattan, plot_manhattan!, manhattancolor
    
   
    # eqtl functions 
    include("./eqtl/eqtl_recipe.jl");
    export eqtlplot, EQTLPlot   

    include("./eqtl/eqtl_helpers.jl");
    export plot_eQTL



end # module BigRiverQTLPlots
