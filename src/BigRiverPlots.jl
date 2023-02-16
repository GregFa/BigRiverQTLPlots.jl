module BigRiverPlots

    # dependent packages 
    using DataFrames, Statistics, Helium
    using RecipesBase
    # using CSV, DelimitedFiles, Missings 
    # using LinearAlgebra, Statistics, Optim
    # using Random, Distributions
    
    # utils functions
    include("./utils.jl");
    export sortnatural, pseudoticks
    
    include("./utils_qtl.jl");
    export plotQTL
    
    include("./utils_eqtl.jl");
    export ploteQTL
    
    # plotting recipes functions 
    include("./recipes/eqtl_recipe.jl");
    export eqtlplot   
    include("./recipes/qtl_recipe.jl");
    export qtlplot   


end # module BigRiverPlots
