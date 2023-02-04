module BigriverPlots

    # dependent packages 
    using Helium
    using RecipesBase
    # using CSV, DelimitedFiles, DataFrames, Missings 
    # using LinearAlgebra, Statistics, Optim
    # using Random, Distributions, LoopVectorization
    
    # utils functions
    include("./utils.jl");
    export sortnatural, pseudoticks
    
    # recipes for plotting 
    # include("./recipes/eqtl_recipe.jl");
    # export eqtlplot   
    include("./recipes/qtl_recipe.jl");
    export qtlplot   


end # module BigriverPlots


