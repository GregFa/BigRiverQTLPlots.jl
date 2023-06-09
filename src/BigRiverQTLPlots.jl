module BigRiverQTLPlots

    # dependent packages 
    using DataFrames, Statistics, Helium
    using RecipesBase
    # using CSV, DelimitedFiles, Missings 
    # using LinearAlgebra, Statistics, Optim
    # using Random, Distributions
    
    # utils functions
    include("utils.jl");
    export sortnatural, pseudoticks,calcKinship2
    

    # qtl functions
    include("./qtl/qtl_recipe.jl");
    export qtlplot, QTLPlot  

    include("./qtl/qtl_helpers.jl");
    export plot_QTL, get_plot_QTL_inputs
    
    
    # eqtl functions 
    include("./eqtl/eqtl_recipe.jl");
    export eqtlplot, EQTLPlot   

    include("./eqtl/eqtl_helpers.jl");
    export plot_eQTL



end # module BigRiverQTLPlots
