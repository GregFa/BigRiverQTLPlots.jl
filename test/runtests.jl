using BigRiverPlots, BulkLMM
using Helium, Statistics
using Plots
using FileIO
using Test


@testset "BigRiverPlots.jl" begin
    include("utils_tests.jl")
    include("recipes_tests.jl")
end