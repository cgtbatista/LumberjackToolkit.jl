#! /bin/bash

SYSTEMS=("dried" "wet" "nolcc")
HEMICELLULOSE=(
    "XY1" "XY2" "XY3" "XY4" "XY5" "XY6" "XY7" "XY8" "XY9" "XY10" "XY11" "XY12" "XY13" "XY14" "XY15" "XY16" "XY21" "XY22" "XA17" "XB17"
    "XC17" "XA18" "XB18" "XC18" "XA19" "XB19" "XC19" "XA20" "XB20" "XC20" "AN1" "AN2" "AN3" "AN4" "AN5" "AN6" "AN7" "AN8" "AN9" "AN10"
    "AN11" "AN12" "AN13" "AN14" "AN15" "AN16" "AN17" "AN18" "AN19" "AN20" "AN21" "AN22" "AN23" "AN24" "AN31" "AN32" "AN33" "MA25" "MB25"
    "MC25" "MA26" "MB26" "MC26" "MA27" "MB27" "MC27" "MA28" "MB28" "MC28" "MA29" "MB29" "MC29" "MA30" "MB30" "MC30"
)

for MODEL in "${SYSTEMS[@]}"
do
    echo "-----------------------------------------------------------------------------------------------------"
    echo "System: $MODEL                                                             $(date)"
    echo "-----------------------------------------------------------------------------------------------------"
    for HC in "${HEMICELLULOSE[@]}"
    do
        echo "Hemicellulose: $HC"
        echo "-----------------------------------------------------------------------------------------------------"
        JULIA_INPUT="t28_${MODEL}_${HC}.jl"
        OUTPUT="t28_${MODEL}_${HC}.dat"

        cat /dev/null > "$JULIA_INPUT"
        {
            echo "using LumberjackToolkit"
            echo ""
            echo "cellulose = ["
            echo "    \"A24\", \"A23\", \"A30\", \"A29\", \"A28\", \"A19\", \"A11\", \"A12\", \"A5\", \"A6\", \"A7\", \"A16\","
            echo "    \"B24\", \"B23\", \"B30\", \"B29\", \"B28\", \"B19\", \"B11\", \"B12\", \"B5\", \"B6\", \"B7\", \"B16\","
            echo "    \"C24\", \"C23\", \"C30\", \"C29\", \"C28\", \"C19\", \"C11\", \"C12\", \"C5\", \"C6\", \"C7\", \"C16\","
            echo "    \"D24\", \"D23\", \"D30\", \"D29\", \"D28\", \"D19\", \"D11\", \"D12\", \"D5\", \"D6\", \"D7\", \"D16\","
            echo "    \"E24\", \"E23\", \"E30\", \"E29\", \"E28\", \"E19\", \"E11\", \"E12\", \"E5\", \"E6\", \"E7\", \"E16\","
            echo "    \"F24\", \"F23\", \"F30\", \"F29\", \"F28\", \"F19\", \"F11\", \"F12\", \"F5\", \"F6\", \"F7\", \"F16\""
            echo "]"
            echo ""
            echo "init_path = \"/home/users/carletoc/PhD/softwoods/unwrap/systems/\""
            echo "pdb = joinpath(init_path, \"$MODEL\", \"$MODEL.pdb\")"
            echo "dcd = joinpath(init_path, \"$MODEL\", \"$MODEL.dcd\")"
            echo "Mcel = mapwater("
            echo "    pdb, dcd,"
            echo "    selection = at -> in(at.segname, cellulose),"
            echo "    water_selection = at -> at.resname == \"TIP3\", cutoff=2.8,"
            echo "    first=1000"
            echo ")"
            echo "Mman = mapwater("
            echo "    pdb, dcd,"
            echo "    selection = at -> at.segname == \"$HC\","
            echo "    water_selection = at -> at.resname == \"TIP3\", cutoff=2.8,"
            echo "    first=1000"
            echo ")"
            echo "M = mapwater(Mcel, Mman, 2)"
            echo "t = t_residence(M, timestep=0.2)"
            echo "open(\"$OUTPUT\", \"w\") do f"
            echo "    for i in t"
            echo "        println(f, \"\$i\")"
            echo "    end"
            echo "end"
        } >> "$JULIA_INPUT"
        julia -t 2 --project=. "$JULIA_INPUT" &
    done
done

wait
echo "-----------------------------------------------------------------------------------------------------"