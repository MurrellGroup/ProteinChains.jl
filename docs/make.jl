using ProteinChains
using Documenter, Literate

LITERATE_INPUT = joinpath(@__DIR__, "..", "examples")
LITERATE_OUTPUT = OUTPUT = joinpath(@__DIR__, "src/generated")

for (root, _, files) ∈ walkdir(LITERATE_INPUT), file ∈ files
    @show root, files
    # ignore non julia files
    splitext(file)[2] == ".jl" || continue
    # full path to a literate script
    ipath = joinpath(root, file)
    # generated output path
    opath = splitdir(replace(ipath, LITERATE_INPUT=>LITERATE_OUTPUT))[1]
    # generate the markdown file calling Literate
    Literate.markdown(ipath, opath)
end

DocMeta.setdocmeta!(
    ProteinChains,
    :DocTestSetup,
    :(using ProteinChains);
    recursive=true,
)

makedocs(;
    modules=[ProteinChains],
    authors="Anton Oresten <anton.oresten42@gmail.com> and contributors",
    sitename="ProteinChains.jl",
    format=Documenter.HTML(;
        canonical="https://MurrellGroup.github.io/ProteinChains.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API" => "api.md",
        "generated/chain.md"
    ],
    doctest=true,
)

deploydocs(;
    repo="github.com/MurrellGroup/ProteinChains.jl",
    devbranch="main",
)
