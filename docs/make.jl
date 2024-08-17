using ProteinChains
using Documenter

DocMeta.setdocmeta!(ProteinChains, :DocTestSetup, :(using ProteinChains); recursive=true)

makedocs(;
    modules=[ProteinChains],
    authors="anton083 <anton.oresten42@gmail.com> and contributors",
    sitename="ProteinChains.jl",
    format=Documenter.HTML(;
        canonical="https://MurrellGroup.github.io/ProteinChains.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/MurrellGroup/ProteinChains.jl",
    devbranch="main",
)
