# ProteinChains

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://MurrellGroup.github.io/ProteinChains.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://MurrellGroup.github.io/ProteinChains.jl/dev/)
[![Build Status](https://github.com/MurrellGroup/ProteinChains.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MurrellGroup/ProteinChains.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/MurrellGroup/ProteinChains.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/MurrellGroup/ProteinChains.jl)

This Julia package provides implements the `ProteinChain` type: an efficient and GPU-friendly representation of protein chains. 

## Installation

```julia
using Pkg
Pkg.add("ProteinChains")
```

## Examples

```julia
julia> using ProteinChains

julia> structure = pdb"1EYE" # convenient macro to download proteins from the PDB
[ Info: Downloading file from PDB: 1EYE
1-chain ProteinStructure "1EYE.cif":
  2 fields:
    name::String = "1EYE.cif"
    chains::Vector{ProteinChain{Float64}} = <exceeds max length>
  2 properties:
    ids::Vector{String} = ["A"]
    lengths::Vector{Int64} = [253]

julia> chain = structure["A"]
253-residue ProteinChain "A":
  4 fields:
    id::String = "A"
    sequence::String = <exceeds max length>
    backbone::Array{Float64,3} = <exceeds max length>
    atoms::Vector{Vector{ProteinChains.Atom{Float64}}} = <exceeds max length>
  2 properties:
    numbering::Vector{Int64} = <exceeds max length>
    modelnum::Int64 = 1

julia> chain.numbering
253-element Vector{Int64}:
   5
   6
   7
   8
   ⋮
 271
 272
 273
 274
```

## See also
- [Backboner.jl](https://github.com/MurrellGroup/Backboner.jl)
- [BioStructures.jl](https://github.com/BioJulia/BioStructures.jl)
- [PDBTools.jl](https://github.com/m3g/PDBTools.jl)
