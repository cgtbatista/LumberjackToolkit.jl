function center_protein(
                    psfname::String,
                    pdbname::String,
                    trajectory::String;
                    selection="protein", pbc_dimensions=nothing, new_trajectory=nothing,
                    reference=0, vmd="vmd", DebugVMD=false
                )
    
    new_trajectory = isnothing(new_trajectory) ? tempname() * ".dcd" : new_trajectory
    
    if isnothing(pbc_dimensions)
        xscname = replace(trajectory, ".dcd" => ".xsc")
        pbc_dimensions = pbc(xscname)
    elseif typeof(pbc_dimensions) == String
        pbc_dimensions = pbc_dimensions
    else
        throw(ArgumentError("The PBC variable must be a string."))
    end

    tcl = tempname() * ".tcl"

    vmdinput = Base.open(tcl, "w")
    
    Base.write(vmdinput,
        raw"""
        package require pbctools

        mol new""" * " $psfname " * raw"""

        mol addfile""" * " $pdbname " * raw"""

        mol addfile""" * " $trajectory waitfor all" * raw"""

        
        animate goto 0

        """ * "$pbc_dimensions" * raw"""

        pbc wrap -centersel""" * " \"$selection\" -center com -compound residue -all" * raw"""
        

        set ref_frame """ * "$reference" * raw"""

        set ref_selection [atomselect top """ * "\"$selection\" frame \$ref_frame]" * raw"""


        set num_frames [molinfo top get numframes]
        for {set i 0} {$i < $num_frames} {incr i} {
            
            set cur_selection [atomselect top """ * "\"$selection\" frame \$i]" * raw"""


            set mat_trans [measure fit $cur_selection $ref_selection]

            set sel_all [atomselect top "all" frame $i]
            $sel_all move $mat_trans
        }

        animate write dcd""" * " $new_trajectory" * raw"""
        """
        )

    Base.close(vmdinput)

    vmdoutput = split(Base.read(`$vmd -dispdev text -e $tcl`, String), "\n")

    if DebugVMD
        return vmdoutput
    else
        return new_trajectory
    end
end

function pbc(xscname::String)

    xscfile = split(Base.read(xscname, String), "\n")

    pbc_string = ""

    for line in xscfile
        
        if !occursin("#", line) && line != ""
            
            l = split(line)
            
            x = [ parse(Float64, l[2]), parse(Float64, l[3]), parse(Float64, l[4]) ]
            y = [ parse(Float64, l[5]), parse(Float64, l[6]), parse(Float64, l[7]) ]
            z = [ parse(Float64, l[8]), parse(Float64, l[9]), parse(Float64, l[10]) ]
            a, b, c = acosd( dot(x, y) / (norm(x) * norm(y)) ), acosd( dot(x, z) / (norm(x) * norm(z)) ), acosd( dot(y, z) / (norm(y) * norm(z)) )
            
            pbc_string = "pbc set { $(sum(x))) $(sum(y)) $(sum(z)) $a $b $c } -all"

        end
    end

    return pbc_string
end