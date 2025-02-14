## VMD -- solvating cellulose macrofibrils

package require psfgen
topology /tmp/jl_KHy12jKrS3.rtf

##pdbalias residue HOH TIP3
readpsf  /home/carlos/Documents/Sandbox/paul/misc/cellulose.psf
coordpdb /home/carlos/Documents/Sandbox/paul/misc/system.pdb

segment WATA {
      auto none
      pdb  /tmp/jl_UJMTOBafoT.pdb
  }
coordpdb /tmp/jl_UJMTOBafoT.pdb WATA

guesscoord

writepsf /home/carlos/Documents/Sandbox/paul/setup/system.psf
writepdb /home/carlos/Documents/Sandbox/paul/setup/system.pdb

exit
