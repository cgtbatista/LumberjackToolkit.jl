"""
    vmd_get_charges(psf::String, pdb::String; newpdb=nothing, pdb_column="beta", vmd="vmd")

Uses VMD to catch the charges inside the PSF file and transfer it to PDB file on b-factor or occupancy columns.
Notice that you need VMD installed to run this code. Check it on: https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD.

# Arguments
- `psf::String` and `pdb::String`: the name of the raw files.
- `newpdb::String`: the name of the new PDB file with the charges. The default is `nothing` and it will generate a PDB on your temporary files.
- `pdb_column::String`: the column to store the charges. The default is `beta`.
- `vmd::String`: the VMD executable. The default is `vmd`.
"""
function vmd_get_charges(psf::String, pdb::String; newpdb=nothing, pdb_column="beta", vmd="vmd")
    
    if (last(split(psf, ".")) != "psf") || (last(split(pdb, ".")) != "pdb")
        throw(ArgumentError("The input files must be on the PSF and PDB format respectively."))
    end

    if Sys.which(vmd) ≡ nothing
        throw(ArgumentError("The VMD executable was not found. Check it on: https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD"))
    end

    if isnothing(newpdb); newpdb = tempname() * ".pdb"; end

    vmdinput_file = tempname() * ".tcl"
    
    vmdinput = open(vmdinput_file, "w")
    Base.write(vmdinput, "mol new \"$psf\" \n")
    Base.write(vmdinput, "mol addfile \"$pdb\" \n")
    Base.write(vmdinput, "set sel [ atomselect top \"all\" ] \n")
    Base.write(vmdinput, "\$sel set $pdb_column [\$sel get charge] \n")
    Base.write(vmdinput, "\$sel writepdb $newpdb \n")
    Base.write(vmdinput, "\n")
    Base.write(vmdinput, "exit")
    Base.close(vmdinput)
    vmdoutput = split(Base.read(`$vmd -dispdev text -e $vmdinput_file`, String), "\n")

    Base.rm(vmdinput_file)
    
    return vmdoutput
end