"""
    vmd_trajectory_alignment(psfname::String, trjname::String; new_trajectory=nothing, selection="not water", reference=0, vmd="vmd", DebugVMD=false)

Align the frames of a trajectory using the VMD `measure fit` command. The selection can be defined by the user,
and the default is `not water`. The reference frame can be defined by the user, and the default is `0`.

### Arguments
- `psfname::String`: The name of the PSF file.
- `trjname::String`: The name of the trajectory file. It can be in any format that VMD can read.
- `new_trajectory::String=nothing`: The name of the new DCD trajectory file. The default creates a temporary file.
- `selection::String="not water"`: The selection to be used to align the frames.
- `reference=0`: The reference frame to be used to align the frames.
- `vmd="vmd"`: The VMD executable. The default is `vmd`.
- `DebugVMD=false`: If `true`, the output of VMD will be printed.
"""
function align_frames(
    psfname::String, trjname::String; new_trajectory=nothing,
    selection="not water", reference::Int64=0,
    vmd="vmd", DebugVMD=false
)    
    new_trajectory = isnothing(new_trajectory) ? tempname() * ".dcd" : new_trajectory
    tcl = tempname() * ".tcl"
    Base.open(tcl, "w") do file
        println(file, """
        package require pbctools

        mol new     $psfname
        mol addfile $trjname waitfor all
        
        animate goto 0
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
        """)
    end
    vmdoutput = Base.read(`$vmd -dispdev text -e $tcl`, String)
    return DebugVMD ? vmdoutput : new_trajectory
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
    trj = writepdb_trajectory(
        pdbname, trajectory; selection=selection, first=first, last=last, step=step
    )
    mapfractions = tempname() * ".out"
    mdlovofit() do exe
        run(pipeline(`$exe -mapfrac $trj`; stdout=mapfractions))
    end
    data = readdlm(mapfractions, comments=true, comment_char='#')
    range = 1:findlast(<(1), data[:,1])
    fraction, RMSDl, RMSDh, RMSD = data[range,1], data[range,2], data[range,3], data[range,4]
    println("""
    -----------------
    Mapping Fractions
    -----------------
    $trj, $mapfractions
    
    Greatest fraction for which the RMSD-low is smaller than 1: $(round(fraction[findlast(<(1.0), RMSDl)],digits=2))

    fraction: contains the fraction of atoms considered in the alignment.
    RMSDl   : contains the RMSD of the fraction of the structure with the lowest RMSD. The best aligned!
    RMSDh   : contains the RMSD of the fraction not considered for the alignment.
    RMSD    : contains the RMSD of the whole structure.
    """)
    return fraction, RMSDl, RMSDh, RMSD
end


function lovo_fitting(
    pdbname::String, trajectory::String, fraction::AbstractFloat;
    selection="protein", fitting="protein and name CA", alignedPDB=nothing, reference=1,
    first=1, last=nothing, step=1
)    
    atoms = PDBTools.readPDB(pdbname, selection)
    tmp_trajectory = writepdb_trajectory(
        pdbname, trajectory;
        selection=selection, first=first, last=last, step=step
    )
    fitatoms = isnothing(fitting) ? atoms : PDBTools.readPDB(pdbname, fitting)
    idx_map = Dict(atom.index => idx for (idx, atom) in enumerate(fitatoms))
    tempPDB_atoms2consider = tempname() * ".pdb"
    for atom in atoms
        idx = get(idx_map, atom.index, nothing)
        if !isnothing(idx)
            atom.beta = 1.0
        else
            atom.beta = 0.0
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