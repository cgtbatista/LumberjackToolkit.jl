function calcular_diametro_fibrila(data; nbins=100)
    xyz = readdlm(data)
    x, y, z = xyz[:,1], xyz[:,2], xyz[:,3]
    zmin, zmax = extrema(z)
    bin_edges = range(
        zmin,
        zmax,
        length=nbins+1
    )
    diameters = Float64[]
    for i in 1:nbins
        lower = bin_edges[i]
        upper = i < nbins ? bin_edges[i+1] : z_max + eps()
        mask = (z .>= lower) .& (z .< upper)
        slice_x = x[mask]
        slice_y = y[mask]
        length(slice_x) < 2 && continue ## skip if there are less than 2 points
        d = 0.0
        n_pontos = length(slice_x)
        for j in 1:n_pontos
            for k in (j+1):n_pontos
                dist = hypot(slice_x[j] - slice_x[k], slice_y[j] - slice_y[k])
                d = max(d, dist)
            end
        end
        
        push!(diameters, d)
    end
    
    # Gerar histograma
    histogram(diameters, bins=20, xlabel="Diâmetro", ylabel="Frequência",
             title="Distribuição de Diâmetros", legend=false)
    
    return (d_mean = mean(diameters),
            d_media = median(diameters),
            d_sd = std(diameters),
            hist = current())
end