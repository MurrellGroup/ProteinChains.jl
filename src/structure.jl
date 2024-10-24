"""
    ProteinStructure{T} <: AbstractVector{ProteinChain{T}}

## Fields
- `name::String`: Usually just the base name of the original file.
- `atoms::Vector{Atom{T}}`: free atoms from the structure that were not part of any protein chain.
- `chains::Vector{<:ProteinChain{T}}`: a collection of `ProteinChain`s.

## Examples

```jldoctest
julia> structure = pdb"1ASS"
```
"""
struct ProteinStructure{T} <: AbstractVector{ProteinChain{T}}
    name::String
    atoms::Vector{Atom{T}}
    chains::Vector{<:ProteinChain{T}}
end

Base.convert(::Type{ProteinStructure{T}}, structure::ProteinStructure) where T =
    ProteinStructure(structure.name, convert(Vector{Atom{T}}, structure.atoms), convert(Vector{ProteinChain{T}}, structure.chains))

Base.size(structure::ProteinStructure) = (length(structure.chains),)

_chainid_to_index(structure::ProteinStructure, id::AbstractString) = findfirst(c -> c.id == id, structure.chains)

Base.getindex(structure::ProteinStructure, i::Integer) = structure.chains[i]
Base.getindex(structure::ProteinStructure, id::AbstractString) = structure[_chainid_to_index(structure, id)]
Base.getindex(structure::ProteinStructure, is::AbstractVector) = ProteinStructure(structure.name, structure.atoms, map(i -> structure[i], is))

Base.setindex!(structure::ProteinStructure, chain::ProteinChain, i::Integer) = (structure.chains[i] = chain)
Base.setindex!(structure::ProteinStructure, chain::ProteinChain, id::AbstractString) = (structure.chains[_chainid_to_index(structure, id)] = chain)

Base.summary(structure::ProteinStructure) = "$(length(structure))-chain $(typeof(structure)) \"$(structure.name)\""

function Base.show(io::IO, structure::ProteinStructure)
    print(io, "$(typeof(structure))(")
    for fieldname in fieldnames(ProteinStructure)
        show(io, getfield(structure, fieldname))
        fieldname != :chains && print(io, ", ")
    end
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", structure::ProteinStructure)
    print(io, summary(structure))
    for chain in structure
        print(io, "\n ", summary(chain))
    end
end

function map_atoms!(f::Function, structure::ProteinStructure, args...)
    for chain in structure
        map_atoms!(f, chain, args...)
    end
    for i in eachindex(structure.atoms)
        structure.atoms[i] = f(structure.atoms[i], args...)
    end
    return structure
end

addproperties(structure::ProteinStructure, names::Symbol...) =
    ProteinStructure(structure.name, structure.atoms, addproperties.(structure.chains, names...))