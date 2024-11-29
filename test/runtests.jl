using LumberjackToolkit
using Test

### DATASETS

## POPC
popc_pdb = joinpath(@__DIR__, "../test/popc/popc.pdb")
popc_psf = joinpath(@__DIR__, "../test/popc/popc.psf")
popc_dcd = joinpath(@__DIR__, "../test/popc/popc.dcd")
# only water and POPC lipid molecules.
# NPT ensemble (300K and 1.0 atm)

## S2 - Softwood Secondary Cell Wall

### MAIN TESTS

@testset "Checking the densityprofile(...)" begin

    dist1, dens1 = densityprofile(popc_pdb, popc_dcd, profile="number", selection="not water", echo=false)
    dist2, dens2 = densityprofile(popc_pdb, popc_dcd, profile="number", selection="not water", echo=false, resolution=2.0)
    dist3, dens3 = densityprofile(popc_pdb, popc_dcd, profile="charge", selection="not water", echo=false, resolution=2.0)
    dist4, dens4 = densityprofile(popc_pdb, popc_dcd, profile="electron", selection="not water", echo=false, resolution=2.0)
    dist5, dens5 = densityprofile(popc_pdb, popc_dcd, profile="mass", selection="not water", echo=false, resolution=2.0)

    #idx = rand(eachindex(dist1[1]), 2)
    idx = rand(1:10, 2); i, j = idx[1], idx[2]

    # matching the output lengths...
    @test length(dist1[i]) == length(dist1[j]) && length(dens1[i]) == length(dens1[j])
    @test length(dist1[i]) != length(dist2[j]) && length(dens1[i]) != length(dens2[j])
    @test length(dist1) == length(dens1) && length(dist2) == length(dens2)
    @test length(dist1) == length(dist2) && length(dens1) == length(dens2)
    @test length(dist2) == length(dist3) == length(dist4) == length(dist5)

    # output content
    @test typeof(dens1) == typeof(dist1) == Vector{Vector{Float64}}
    @test dens2 != dens3 != dens4 != dens5
    @test dens1 != dens2 && dens1[i] != dens1[j] && dens2[i] != dens2[j]
    @test dist1 != dist2 && dist1[i] != dist1[j] && dist2[i] != dist2[j]

    # Write your tests here.
end
