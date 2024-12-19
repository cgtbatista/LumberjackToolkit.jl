#mol new ../setup/sys.psf
#mol addfile ../all_5dried.dcd

##

function rmsd_pcw(filename::String; tcl_script=nothing)

    if isnothing(tcl_script); tcl_script = tempname() * ".tcl"; end

    vmdinput = open(tcl_script, "w")

    Base.write(vmdinput, "mol new \"$filename.psf\" \n")
    Base.write(vmdinput, "mol addfile \"$filename.pdb\" \n")
    Base.write(vmdinput, "mol addfile \"$filename.dcd\" waitfor -1 \n")
    Base.write(vmdinput, "\n")
    Base.write(vmdinput, "animate goto 0 \n")
    Base.write(vmdinput, "set sel [atomselect top \"not water and not ions\"] \n")



    Base.write(vmdinput, "\n")
    Base.write(vmdinput, "exit")
    Base.close(vmdinput)

    vmdoutput = split(Base.read(`$vmd -dispdev text -e $vmdinput_file`, String), "\n")
    
    return vmdoutput, filename
end


