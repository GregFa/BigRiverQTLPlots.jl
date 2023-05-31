module BigRiverPlots

    # dependent packages 
    using DataFrames, Statistics, Helium
    using RecipesBase
    # using CSV, DelimitedFiles, Missings 
    # using LinearAlgebra, Statistics, Optim
    # using Random, Distributions
    
    # utils functions
    include("./utils/utils.jl");
    export sortnatural, pseudoticks
    
    include("./utils/utils_qtl.jl");
    export plot_QTL, get_plot_QTL_inputs
    
    include("./utils/utils_eqtl.jl");
    export plot_eQTL
    
    # plotting recipes functions 
    include("./recipes/eqtl_recipe.jl");
    export eqtlplot, EQTLPlot   
    include("./recipes/qtl_recipe.jl");
    export qtlplot, QTLPlot   


end # module BigRiverPlots
