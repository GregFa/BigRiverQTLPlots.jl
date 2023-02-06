using BigriverPlots
using Helium
using Test

# Generate test data
include("test_data.jl")

@testset "BigriverPlots.jl" begin
    include("utils_test.jl")
end