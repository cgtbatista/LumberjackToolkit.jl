module LumberjackToolkit

import PDBTools
import EasyFit
import MolSimToolkit
import StatsPlots, Plots
import Interpolations

import FFTW: rfft, irfft
import MDLovoFit_jll: mdlovofit
import DelimitedFiles: readdlm, writedlm
import Printf: @sprintf, @printf
import StaticArrays: @SVector, SVector, SMatrix, @SMatrix, MArray
import Statistics: mean, median, std
import StatsBase: Histogram, fit
import LinearAlgebra: norm, dot, cross, diag, diagm, eigen, det, normalize, BLAS.set_num_threads
import Base.Threads: @threads, nthreads

## distances -- compute the distances between two sets of atoms such as: cellulose-water, cellulose-hemicellulose, catalytic site, etc.
export min_distances

## RMSD
export rmsd_plot

export avg_densityprofile

export CarbohydrateDihedrals
export simulation_steps, dihedral_atoms, dihedral_indexes, dihedrals
export get_pressure, get_energy

## PMF -- Potential of Mean Force
export WHAM, pmf

export pmf1D, pmf2D

export align_frames, namd_pbc
export center_of_mass, rmsd

export write_frame!, writepdb_trajectory, lovo_mapping, lovo_fitting
export pdb2trajectory, frame_positions

export readcoords

## PCA
export PCA, _pca_data
export r_avg, dispm, covarm, massm, mass, pca, projections, explained_variance, quasiharmonic_frequency, quasiharmonic_modes
export E_quasiharmonic, P_quasiharmonic, pmf_pca

##
export simsteps, realtime
export testfiles, diffusion, msd, frame_coordinates, molindexes, writecoords, molframes

## TEST!!
export readcoords2, msd2, unwrap, dimcells, displace, originaldisplace
export teste2, diameter_analysis, fibrilradii, chain_centers

export mapwater, checking_residence, t_residence, closest2fibril, water_hbonding, residence, water_hbonding_parallel

## dupree
export cellulose_surface, filterSTL, fibrilwidth, fibril_surface, fibril_slice, catalytic_distances
export chargesPSF, binning, binspecs, align_frames, electrons, density_profile, notaxis, axis2dims

# Getting some properties profiles using the density distribution (e.g. electron density profile of POPC on the box)
include("./residence_time.jl")

include("./density.jl")

include("./dihedrals.jl")

include("./plotting.jl")

include("./utils.jl")

include("./PMF.jl")

include("./align.jl")

include("./misc-namd.jl")
include("./pca.jl")

#
include("./diffusion.jl")
include("./simulation.jl")

include("./unwrap.jl")

include("./dupree.jl")
include("./distances.jl")

end