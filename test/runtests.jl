using LumberjackToolkit
using Test

### MAIN TESTS
@testset "density.jl -- TRUE POSITIVES" begin
    # Test 1: Matching the vector lengths
    dist, dens = densityprofile("popc/popc.pdb", "popc/popc.dcd", profile="number", selection="not water")
    @test length(dist) == length(dens)
    # Write your tests here.
end
