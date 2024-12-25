"""
    align_trajectory(psfname::String, pdbname::String, trajectory::String; ...)

Aligns a trajectory to a reference structure using VMD. It loads the psf, pdb, and trajectory files, wraps the trajectory, and aligns the frames to the reference structure.
    The `selection` argument is used to define the atoms to be aligned based on the VMD selection, such as a protein or a cellulose fibril.
    The `DebugVMD` argument equals  **true** will return  the VMD output. The function returns the path to the new trajectory.
    If the `new_trajectory` argument is not provided, a temporary file is created.

# Arguments
- `psfname::String`: The path to the PSF file.
- `pdbname::String`: The path to the PDB file.
- `trajectory::String`: The path to the trajectory file.
- `selection::String`: The VMD selection to be aligned. Default is "protein".
- `pbc_dimensions::String`: The PBC dimensions to be set on **pbctools format**. If not provided, the function will try to read the dimensions from the .xsc file, based on the location of trajectory.
- `new_trajectory::String`: The path to the new trajectory file. If not provided, a temporary file is created.
- `reference::Int`: The frame to be used as reference. Default is 0, the initial frame.
- `vmd::String`: The path to the VMD executable. Default is "vmd".
- `DebugVMD::Bool`: If **true**, the function will return the VMD output. Default is **false**.
"""
function align_trajectory(
                    psfname::String,
                    pdbname::String,
                    trajectory::String;
                    selection="protein and name CA", pbc_dimensions=nothing, new_trajectory=nothing, reference=0, vmd="vmd", DebugVMD=false
                )
    
    new_trajectory = isnothing(new_trajectory) ? tempname() * ".dcd" : new_trajectory
    
    if isnothing(pbc_dimensions)
        xscname = replace(trajectory, ".dcd" => ".xsc")
        pbc_dimensions = namd_pbc(xscname)
    elseif typeof(pbc_dimensions) == String
        pbc_dimensions = pbc_dimensions
    else
        throw(ArgumentError("The PBC variable must be a string."))
    end

    tcl = tempname() * ".tcl"

    vmdinput = Base.open(tcl, "w")
    
    Base.write(vmdinput,"""
        package require pbctools

        mol new     $psfname
        mol addfile $pdbname
        mol addfile $trajectory waitfor all
        
        animate goto 0
        $pbc_dimensions
        pbc wrap -centersel \"$selection\" -center com -compound residue -all
        
        set ref_frame $reference
        set ref_selection [atomselect top \"$selection\" frame \$ref_frame]

        set num_frames [molinfo top get numframes]

        for {set i 0} {\$i < \$num_frames} {incr i} {
            set cur_selection [atomselect top \"$selection\" frame \$i]
            set mat_trans [measure fit \$cur_selection \$ref_selection]
            set sel_all [atomselect top "all" frame \$i]
            \$sel_all move \$mat_trans
        }

        animate write dcd $new_trajectory
        """
        )

    Base.close(vmdinput)

    vmdoutput = split(Base.read(`$vmd -dispdev text -e $tcl`, String), "\n")

    if DebugVMD
        return vmdoutput, tcl
    else
        return new_trajectory
    end
end

"""
    map_fractions(atoms::AbstractVector{<:PDBTools.Atom}, trajectory_file::String)


"""
function lovo_mapping(
                    pdbname::String, 
                    trajectory::String;
                    selection="protein and name CA",
                    first=1, last=nothing, step=1
                )

    tmp_trajectory = writepdb_trajectory(pdbname, trajectory; selection=selection, first=first, last=last, step=step)

    mapfractions = tempname() * ".out"

    mdlovofit() do exe
        run(pipeline(`$exe -mapfrac $tmp_trajectory`; stdout=mapfractions))
    end

    data = readdlm(mapfractions, comments=true, comment_char='#')
    range = 1:findlast(<(1), data[:,1])
    fraction, RMSDl, RMSDh, RMSD = data[range,1], data[range,2], data[range,3], data[range,4]

    println("""
    -----------------
    Mapping Fractions
    -----------------
    $tmp_trajectory, $mapfractions
    
    Greatest fraction for which the RMSD-low is smaller than 1: $(round(fraction[findlast(<(1.0), RMSDl)],digits=2))

    fraction: contains the fraction of atoms considered in the alignment.
    RMSDl   : contains the RMSD of the fraction of the structure with the lowest RMSD. The best aligned!
    RMSDh   : contains the RMSD of the fraction not considered for the alignment.
    RMSD    : contains the RMSD of the whole structure.
    """)

    return fraction, RMSDl, RMSDh, RMSD
end



function lovo_fitting(
                    pdbname::String, 
                    trajectory::String,
                    fraction::AbstractFloat;
                    selection="protein and name CA", alignedPDB=nothing, atoms2consider=nothing, reference=1,
                    first=1, last=nothing, step=1
                )
    
    atoms = PDBTools.readPDB(pdbname, selection)
    tmp_trajectory = writepdb_trajectory(pdbname, trajectory; selection=selection, first=first, last=last, step=step)
    
    if isnothing(atoms2consider)
        atoms2consider = atoms
    elseif !isnothing(atoms2consider) && (typeof(atoms2consider) != AbstractVector{<:PDBTools.Atom})
        throw(ArgumentError("The atoms2consider variable must be an AbstractVector{<:PDBTools.Atom}."))
    end

    idx_map = Dict(atom.index => idx for (idx, atom) in enumerate(atoms2consider))

    tempPDB_atoms2consider = tempname() * ".pdb"
    for atom in atoms
        idx = get(idx_map, atom.index, nothing)
        if !isnothing(idx)
            atoms[idx].beta = 1.0
        else
            atoms[idx].beta = 0.0
        end
    end
    PDBTools.writePDB(atoms, tempPDB_atoms2consider)

    rmsf_output, fitting_output = tempname() * ".rmsf", tempname() * ".out"
    
    alignedPDB = isnothing(alignedPDB) ? tempname() * ".pdb" : alignedPDB
    
    try
        mdlovofit() do exe
            run(pipeline(`$exe -f $fraction -iref $reference -rmsf $rmsf_output -t $alignedPDB $tmp_trajectory`; stdout=fitting_output))
        end
    catch 
        "ERROR in MDLovoFit execution"
        "Command executed: $command"
    end

    data = readdlm(fitting_output; comments=true, comment_char='#')
    
    RMSDl, RMSDh, RMSD = data[:,2], data[:,3], data[:,4]
    RMSF = readdlm(rmsf_output; comments=true, comment_char='#')[:,2]

    println("""
    ------------
    LOVO Fitting
    ------------
    $alignedPDB

    Average RMSD of all atoms: $(round((mean(RMSD)), digits=2))
    Average RMSD of the $(round(100*fraction,digits=1))% atoms of lowest RMSD: $(round((mean(RMSDl)), digits=2))
    Average RMSD of the $(round(100*(1-fraction),digits=1))% atoms of highest RMSD: $(round((mean(RMSDh)), digits=2))

    RMSF data is based on $(length(RMSF)) atoms.
    """)

    return RMSDl, RMSDh, RMSD, RMSF
end