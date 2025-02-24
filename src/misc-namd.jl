function namd_pbc(xscname::String)
    xscfile = split(Base.read(xscname, String), "\n")
    for line in xscfile
        if !occursin("#", line) && line != ""
            l = split(line)         
            x = [ parse(Float64, l[2]), parse(Float64, l[3]), parse(Float64, l[4]) ]
            y = [ parse(Float64, l[5]), parse(Float64, l[6]), parse(Float64, l[7]) ]
            z = [ parse(Float64, l[8]), parse(Float64, l[9]), parse(Float64, l[10]) ]
            a = acosd( dot(x, y) / (norm(x) * norm(y)) )
            b = acosd( dot(x, z) / (norm(x) * norm(z)) )
            c = acosd( dot(y, z) / (norm(y) * norm(z)) )
            return "pbc set { $(sum(x)) $(sum(y)) $(sum(z)) $a $b $c } -all"
        end
    end
end

function get_pressure(logname::String; flag="PRESSURE", unit="bar", last=true)
    return get_pressure(split(Base.read(logname, String), "\n"); flag=flag, unit=unit, last=last)
end

function get_pressure(logfile::Vector{SubString{String}}; flag="PRESSURE", unit="bar", last=true)

    p = Vector{SMatrix}()

    for line in logfile
        l = split(line)

        if !occursin("Info:", line) && startswith(line, flag) && length(l) == 11
            push!(p, @SMatrix(
                [
                    parse(Float64, l[3]) parse(Float64, l[4]) parse(Float64, l[5]);
                    parse(Float64, l[4]) parse(Float64, l[9]) parse(Float64, l[10]);
                    parse(Float64, l[5]) parse(Float64, l[10]) parse(Float64, l[11])
                ])
                )
        end
    end

    if p == []
        throw(ArgumentError("The pressure flag was not recognized. Try `PRESSURE, `GPRESSURE`, `PRESSAVG` or `GPRESSAVG`."))
    end

    if last
        return p[end]
    else
        return p
    end

end


function get_energy(logname::String; property="TEMP", last=true)
    return get_energy(split(Base.read(logname, String), "\n"); property=property, last=last)
end

function get_energy(logfile::Vector{SubString{String}}; property="TEMP", last=true)

    E = Vector{Float64}()

    idx = 13
    idx = lowercase(property) == "bond" ? 3 : idx
    idx = lowercase(property) == "angle" ? 4 : idx
    idx = lowercase(property) == "dihed" ? 5 : idx
    idx = lowercase(property) == "imprp" ? 6 : idx
    idx = lowercase(property) == "elect" ? 7 : idx
    idx = lowercase(property) == "vdw" ? 8 : idx
    idx = lowercase(property) == "boundary" ? 9 : idx
    idx = lowercase(property) == "misc" ? 10 : idx
    idx = lowercase(property) == "kinetic" ? 11 : idx
    idx = lowercase(property) == "total" ? 12 : idx
    idx = lowercase(property) == "temp" ? 13 : idx
    idx = lowercase(property) == "potential" ? 14 : idx
    idx = lowercase(property) == "totalavg" ? 15 : idx
    idx = lowercase(property) == "tempavg" ? 16 : idx
    idx = lowercase(property) == "pressure" ? 17 : idx
    idx = lowercase(property) == "gpressure" ? 18 : idx
    idx = lowercase(property) == "volume" ? 19 : idx
    idx = lowercase(property) == "pressavg" ? 20 : idx
    idx = lowercase(property) == "gpressavg" ? 21 : idx

    for line in logfile
        l = split(line)

        if !occursin("Info:", line) && startswith(line, "ENERGY:")
            push!(E, parse(Float64, l[idx]))
        end
    end

    if E == []
        throw(ArgumentError("The energy flag was not recognized. Try..."))
    end

    if last
        return E[end]
    else
        return E[begin+1:end]
    end

end