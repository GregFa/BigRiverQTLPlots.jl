using BigRiverQTLPlots, BulkLMM
using Statistics, Random
using Plots
using FileIO, Helium
using Test
ENV["GKSwstype"] = "100"

@testset "BigRiverQTLPlots.jl" begin
    include("utils_tests.jl")
    include("recipes_tests.jl")
end