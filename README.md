# ProteinChains

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://MurrellGroup.github.io/ProteinChains.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://MurrellGroup.github.io/ProteinChains.jl/dev/)
[![Build Status](https://github.com/MurrellGroup/ProteinChains.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MurrellGroup/ProteinChains.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/MurrellGroup/ProteinChains.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/MurrellGroup/ProteinChains.jl)

This Julia package provides implements the `ProteinChain` type: a GPU-friendly structure-of-arrays representation of protein chains.

## Installation

```julia
using Pkg
Pkg.add("ProteinChains")
```

## Examples

The `ProteinChain` type is meant to only store a set of quintessential fields, from which most other properties can be derived.

```julia
julia> using ProteinChains

julia> structure = pdb"1EYE" # string macro to fetch proteins from the PDB
[ Info: Downloading file from PDB: 1EYE
1-chain ProteinStructure "1EYE.cif"
 256-residue ProteinChain{Float64} (A)

julia> propertynames(chain)
(:id, :sequence, :backbone, :numbering, :atoms)
```

To store additional properties, `AnnotatedProteinChain` can be used to add dynamic properties to the chain:

```julia
julia> annotated_chain = annotate(chain; model=1)
256-residue AnnotatedProteinChain{Float64} (A):
  6 fields:
    id::String = "A"
    sequence::String = <exceeds max length>
    backbone::Array{Float64,3} = <exceeds max length>
    numbering::Vector{Int64} = <exceeds max length>
    atoms::Vector{Vector{ProteinChains.Atom{Float64}}} = <exceeds max length>
    indexable_properties::Vector{Symbol} = Symbol[]
  1 property:
    model::Int64 = 1
```

For properties of type `<:AbstractArray` that represent residue-level information, `annotate_indexable!` will index the last dimension of the property when the chain is indexed:

```julia
julia> annotate_indexable!(annotated_chain; secondary_structure=assign_secondary_structure(annotated_chain)
256-residue AnnotatedProteinChain{Float64} (A):
  6 fields:
    id::String = "A"
    sequence::String = <exceeds max length>
    backbone::Array{Float64,3} = <exceeds max length>
    numbering::Vector{Int64} = <exceeds max length>
    atoms::Vector{Vector{ProteinChains.Atom{Float64}}} = <exceeds max length>
    indexable_properties::Vector{Symbol} = [:secondary_structure]
  2 properties:
    model::Int64 = 1
    secondary_structure::Vector{Int64} = [1, 1, 3, 3, 3, 3, 3, 3, 3, 1  …  2, 2, 2, 2, 2, 2, 2, 1, 1, 1]
```

## See also

- [Backboner.jl](https://github.com/MurrellGroup/Backboner.jl)
- [BioStructures.jl](https://github.com/BioJulia/BioStructures.jl)
- [PDBTools.jl](https://github.com/m3g/PDBTools.jl)
