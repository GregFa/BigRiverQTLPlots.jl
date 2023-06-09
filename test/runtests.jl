using BigRiverQTLPlots, BulkLMM
using Statistics, Random
using Plots
using FileIO, Helium
using Test


@testset "BigRiverQTLPlots.jl" begin
    include("utils_tests.jl")
    include("recipes_tests.jl")
end