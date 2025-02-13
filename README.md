# ProteinChains

[![Dev](https://img.shields.io/badge/docs-stable-blue.svg)](https://MurrellGroup.github.io/ProteinChains.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://MurrellGroup.github.io/ProteinChains.jl/dev/)
[![Build Status](https://github.com/MurrellGroup/ProteinChains.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MurrellGroup/ProteinChains.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/MurrellGroup/ProteinChains.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/MurrellGroup/ProteinChains.jl)

This Julia package provides implements the `ProteinChain` type: a chain-level structure-of-arrays type representation of proteins, with support for indexing by residue index.

## Installation

```julia
using Pkg
Pkg.add("ProteinChains")
```

## Quickstart

The `ProteinChain` type is meant to only store a basic set of fields, from which some other properties might be derived.

```julia
julia> using ProteinChains

julia> structure = pdb"1EYE" # string macro to fetch proteins from the PDB
[ Info: File exists: 1EYE
1-element ProteinStructure{Float64} "1EYE.cif"
 256-residue ProteinChain{Float64} (A)

julia> chain = structure["A"]
256-residue ProteinChain{Float64} (A)

julia> propertynames(chain)
(:id, :atoms, :sequence, :numbering, :ins_codes)
```

To store additional properties, simply set them like you would set a field of a mutable struct:

```julia
julia> chain.taxid = 83332
83332
```

For chain-wise properties, the values can be wrapped with `Indexable` to index the last dimension of the array when the chain gets indexed:

```julia
julia> chain.rand3 = Indexable(rand(3,256)) # last dimension matches chain length
Indexable wrapping a Matrix{Float64} (get value with `unwrap(x)`)

julia> chain[1:100].rand3 |> unwrap
3×100 view(::Matrix{Float64}, :, 1:100) with eltype Float64:
 0.273545  0.639173  0.92708   …  0.459441  0.196407  0.880034       
 0.981498  0.70263   0.279264     0.552049  0.89274   0.0328866      
 0.169268  0.117848  0.732741     0.301921  0.187094  0.281187

julia> propertynames(chain)
(:id, :atoms, :sequence, :numbering, :ins_codes, :taxid, :rand3)
```

## See also

- [BioStructures.jl](https://github.com/BioJulia/BioStructures.jl)
- [Backboner.jl](https://github.com/MurrellGroup/Backboner.jl)
- [PDBTools.jl](https://github.com/m3g/PDBTools.jl)
