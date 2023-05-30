using BigRiverPlots, BulkLMM
using Helium, Statistics
using Plots
using FileIO
using Test

# Generate test data
# include("test_data.jl")

@testset "BigRiverPlots.jl" begin
    include("utils_tests.jl")
    include("recipes_tests.jl")
end