using IntegratedWienerProcesses
using Documenter

DocMeta.setdocmeta!(IntegratedWienerProcesses, :DocTestSetup, :(using IntegratedWienerProcesses); recursive=true)

makedocs(;
    modules=[IntegratedWienerProcesses],
    authors="Filip Tronarp <filip.tronarp@uni-tuebingen.de> and contributors",
    repo="https://github.com/filtron/IntegratedWienerProcesses.jl/blob/{commit}{path}#{line}",
    sitename="IntegratedWienerProcesses.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://filtron.github.io/IntegratedWienerProcesses.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/filtron/IntegratedWienerProcesses.jl",
    devbranch="main",
)
