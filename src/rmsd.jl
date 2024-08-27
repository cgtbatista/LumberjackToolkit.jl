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


# ## Completo
# set sel [atomselect top "all"]
# # rmsd calculation loop
# set outfile [open all_rmsd.dat w]
# set nf [molinfo top get numframes]
# set frame0 [atomselect top "all" frame 0]
# for {set i 1 } {$i < $nf } { incr i } {
#     $sel frame $i
#     $sel move [measure fit $sel $frame0]
#     puts "RMSD: $i [measure rmsd $sel $frame0]"
# }
# close $outfile

# ### Celulose
# set fibril_a "A11 A19 A12 A20 A28 A5 A13 A21 A29 A6 A14 A22 A30 A7 A15 A23 A16 A24"
# set fibril_b "B11 B19 B12 B20 B28 B5 B13 B21 B29 B6 B14 B22 B30 B7 B15 B23 B16 B24"
# set fibril_c "C11 C19 C12 C20 C28 C5 C13 C21 C29 C6 C14 C22 C30 C7 C15 C23 C16 C24"
# set fibril_d "D11 D19 D12 D20 D28 D5 D13 D21 D29 D6 D14 D22 D30 D7 D15 D23 D16 D24"
# set fibril_e "E22 E30 E7 E15 E23 E16 E24 E11 E19 E12 E20 E28 E5 E13 E21 E29 E6 E14"
# set fibril_f "F11 F19 F12 F20 F28 F5 F13 F21 F29 F6 F14 F22 F30 F7 F15 F23 F16 F24"
# set fibril_g "G11 G19 G12 G20 G28 G5 G13 G21 G29 G6 G14 G22 G30 G7 G15 G23 G16 G24"
# set outfile [open cel_rmsd.dat w];                                             
# set nf [molinfo top get numframes]
# set frame0 [atomselect top "segid $fibril_a $fibril_b $fibril_c $fibril_d $fibril_e $fibril_f $fibril_g" frame 0]
# set sel [atomselect top "segid $fibril_a $fibril_b $fibril_c $fibril_d $fibril_e $fibril_f $fibril_g"]
# # rmsd calculation loop
# for {set i 1 } {$i < $nf } { incr i } {
#     $sel frame $i
#     $sel move [measure fit $sel $frame0]
#     puts $outfile "[measure rmsd $sel $frame0]"
# }
# close $outfile

# ## Xilanos
# set close_xyl "XY1 XY2 XY3 XY4 XY5 XY6 XY7 XY8 XY9 XY10 XY11 XY12 XY13 XY14 XY15 XY16"
# set far_xyl "XY21 XY22 XA17 XB17 XC17 XA18 XB18 XC18 XA19 XB19 XC19 XA20 XB20 XC20"
# set outfile [open xyl_rmsd.dat w];                                             
# set nf [molinfo top get numframes]
# set frame0 [atomselect top "segid $close_xyl $far_xyl" frame 0]
# set sel [atomselect top "segid $close_xyl $far_xyl"]
# # rmsd calculation loop
# for {set i 1 } {$i < $nf } { incr i } {
#     $sel frame $i
#     $sel move [measure fit $sel $frame0]
#     puts $outfile "[measure rmsd $sel $frame0]"
# }
# close $outfile

# ## Mananos
# set close_man "AN1 AN2 AN3 AN4 AN5 AN6 AN7 AN8 AN9 AN10 AN11 AN12 AN13 AN14 AN15 AN16 AN17 AN18 AN19 AN20 AN21 AN22 AN23 AN24 AN31 AN32"
# set far_man "AN33 MA25 MB25 MC25 MA26 MB26 MC26 MA27 MB27 MC27 MA28 MB28 MC28 MA29 MB29 MC29 MA30 MB30 MC30"
# set outfile [open man_rmsd.dat w]
# set nf [molinfo top get numframes]
# set frame0 [atomselect top "segid $close_man $far_man" frame 0]
# set sel [atomselect top "segid $close_man $far_man"]
# # rmsd calculation loop
# for {set i 1 } {$i < $nf } { incr i } {
#     $sel frame $i
#     $sel move [measure fit $sel $frame0]
#     puts $outfile "[measure rmsd $sel $frame0]"
# }
# close $outfile

# ## Ligninas
# set lig_1 "L1 L3 L6 L7 L8 L9 L10 L13 L17 L20 L26 L28 L34 L42 L49 L51 L52 L53 L57 L65 L69 L75 L77 L79 L82 L83 L88 L90 L95 L96 L21 L24 L30 L31"
# set lig_2 "L35 L36 L38 L39 L41 L43 L47 L54 L61 L63 L72 L80 L81 L12A L12B L12C L15A L15B L15C L19A L19B L19C L23A L23B L23C L25A L25B L25C L27A"
# set lig_3 "L27B L27C L29A L29B L29C L32A L32B L32C L44A L44B L44C L45A L45B L45C L46A L46B L46C L56A L56B L56C L64A L64B L64C L66A L66B L66C L68A"
# set lig_4 "L68B L68C L70A L70B L70C L85A L85B L85C L87A L87B L87C L92A L92B L92C L98A L98B L98C K1 K3 K6 K7 K8 K9 K10 K13 K17 K20 K26 K28 K34 K42"
# set lig_5 "K49 K51 K52 K53 K57 K65 K69 K75 K77 K78 K79 K82 K83 K90 K95 K96 K4 K16 K21 K24 K30 K31 K35 K36 K37 K33 K38 K39 K41 K43 K47 K50 K54 K58"
# set lig_6 "K59 K61 K62 K63 K67 K72 K80 K81 K86 K93 K94 K97 J1 J3 J6 J7 J8 J9 J10 J13 J17 J20 J26 J28 J34 J42 J49 J51 J52 J53 J57 J65 J69 J75 J77 J78"
# set lig_7 "J79 J82 J83 J90 J95 J96 J35 J36 J37 J33 J38 J41 J43 J47 J50 J54 J58 J61 J63 J72 J80 J81 J86 J93 J94 J97"
# set outfile [open lig_rmsd.dat w];                                             
# set nf [molinfo top get numframes]
# set frame0 [atomselect top "segid $lig_1 $lig_2 $lig_3 $lig_4 $lig_5 $lig_6 $lig_7" frame 0]
# set sel [atomselect top "segid $lig_1 $lig_2 $lig_3 $lig_4 $lig_5 $lig_6 $lig_7"]
# # rmsd calculation loop
# for {set i 1 } {$i < $nf } { incr i } {
#     $sel frame $i
#     $sel move [measure fit $sel $frame0]
#     puts $outfile "[measure rmsd $sel $frame0]"
# }
# close $outfile

# exit