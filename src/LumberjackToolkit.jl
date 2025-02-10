module LumberjackToolkit

import PDBTools
import EasyFit
import MolSimToolkit
import StatsPlots, Plots

import FFTW: rfft, irfft
import MDLovoFit_jll: mdlovofit
import DelimitedFiles: readdlm
import Printf: @sprintf, @printf
import StaticArrays: @SVector, SVector, SMatrix, @SMatrix, MArray
import Statistics: mean, median, std
import LinearAlgebra: norm, dot, cross, diag, diagm, eigen, det

import Base.Threads: @threads
## RMSD
export rmsd_plot

## Density profiles
export densityprofile
export ρ, δ
export averaging_profile, average_bins, binning, _ordering_symmetry, _correct_center, _correct_symmetry, _get_reference
export vmd_get_charges

export avg_densityprofile

export CarbohydrateDihedrals
export simulation_steps, dihedral_atoms, dihedral_indexes, dihedrals
export get_pressure, get_energy

export WHAM

export pmf1D, pmf2D

export align_trajectory, namd_pbc
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
export testfiles, coef_diffusion, msd, frame_coordinates, molindexes, writecoords, molframes

## TEST!!
export readcoords2, msd2, unwrap, dimcells, displace, originaldisplace
export teste2, diameter_analysis, fibrilradii, chain_centers

export mapwater, checking_residence, t_residence

## dupree
export cellulose_surface, filterSTL, fibrilwidth, fibril_surface, fibril_slice

# Getting some properties profiles using the density distribution (e.g. electron density profile of POPC on the box)
include("./residence_time.jl")
include("./density.jl")
include("./dihedrals.jl")

include("./plotting.jl")
# others
include("./vmd-scripts.jl")

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

end