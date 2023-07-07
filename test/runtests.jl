using BigRiverQTLPlots, BulkLMM
using Plots
using FileIO, Helium
using Test
ENV["GKSwstype"] = "nul"

@testset "BigRiverQTLPlots.jl" begin
    include("utils_tests.jl")
    include("recipes_tests.jl")
end