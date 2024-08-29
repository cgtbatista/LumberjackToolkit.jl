using LumberjackToolkit
using Test, Testing

### MAIN TESTS
@testset "density.jl -- TRUE POSITIVES" begin
    # Test 1: Matching the vector lengths
    dist, dens = densityprofile(popc_pdb, popc_pdb, profile="number", selection="not water")
    @test length(dist) == length(dens)
    # Write your tests here.
end
