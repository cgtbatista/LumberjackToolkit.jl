using LumberjackToolkit
using Documenter

DocMeta.setdocmeta!(LumberjackToolkit, :DocTestSetup, :(using LumberjackToolkit); recursive=true)

makedocs(;
    modules=[LumberjackToolkit],
    authors="cgtbatista <c203748@dac.unicamp.br> and contributors",
    sitename="LumberjackToolkit.jl",
    format=Documenter.HTML(;
        canonical="https://cgtbatista.github.io/LumberjackToolkit.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/cgtbatista/LumberjackToolkit.jl",
    devbranch="main",
)
