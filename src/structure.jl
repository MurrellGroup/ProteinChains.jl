@dynamic mutable struct ProteinStructure{T} <: AbstractVector{ProteinChain{T}}
    name::String
    atoms::Vector{Atom{T}}
    chains::Vector{ProteinChain{T}}
end

"""
    ProteinStructure{T} <: AbstractVector{ProteinChain{T}}

## Fields
- `name::String`: Usually just the base name of the original file.
- `atoms::Vector{Atom{T}}`: free atoms from the structure that were not part of any protein residue.
- `chains::Vector{ProteinChain{T}}`: a collection of `ProteinChain`s.
"""
ProteinStructure

"""
    ProteinStructure(name, chains; properties...)
"""
ProteinStructure(name::AbstractString, chains::Vector{ProteinChain{T}}, args...; kwargs...) where T =
    ProteinStructure(name, Atom{T}[], chains, args...; kwargs...)

Base.convert(::Type{ProteinStructure{T}}, structure::ProteinStructure) where T =
    ProteinStructure(
        structure.name,
        convert(Vector{Atom{T}}, structure.atoms),
        convert(Vector{ProteinChain{T}}, structure.chains);
        propertypairs(structure, NoFields())...,
    )

Base.size(structure::ProteinStructure) = size(structure.chains)

function chainid_to_index(structure::ProteinStructure, id::AbstractString)
    index = findfirst(c -> c.id == id, structure.chains)
    index === nothing && throw(KeyError(id))
    return index
end

Base.getindex(structure::ProteinStructure, i::Integer) = structure.chains[i]
Base.getindex(structure::ProteinStructure, id::AbstractString) = structure[chainid_to_index(structure, id)]

function Base.getindex(structure::ProteinStructure, i)
    args = (structure.name, structure.atoms, structure.chains[i])
    kwargs = Iterators.map(propertypairs(structure, NoFields())) do (name, value)
        name => if value isa AbstractProperty
            checkproperty(structure, value)
            value isa Indexable ? value[i] : value
        else
            value
        end
    end
    return ProteinStructure(args...; kwargs...)
end

Base.setindex!(structure::ProteinStructure, chain::ProteinChain, i::Integer) = (structure.chains[i] = chain)
Base.setindex!(structure::ProteinStructure, chain::ProteinChain, id::AbstractString) = (structure.chains[chainid_to_index(structure, id)] = chain)

Base.showarg(io::IO, structure::ProteinStructure, ::Bool) = print(io, "$(ProteinStructure) \"$(structure.name)\"")

function map_chains!(f::Function, structure::ProteinStructure, args...)
    for (i, chain) in enumerate(structure)
        structure[i] = f(chain, args...)
    end
    return structure
end

function map_atoms!(f::Function, structure::ProteinStructure, args...)
    map_chains!(c -> map_atoms!(f, c, args...), structure)
    for i in eachindex(structure.atoms)
        structure.atoms[i] = f(structure.atoms[i], args...)
    end
    return structure
end
