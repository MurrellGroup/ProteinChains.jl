"""
    ProteinStructure{T} <: AbstractVector{ProteinChain{T}}

## Fields
- `name::String`: Usually just the base name of the original file.
- `atoms::Vector{Atom{T}}`: free atoms from the structure that were not part of any protein residue.
- `chains::Vector{ProteinChain{T}}`: a collection of `ProteinChain`s.
- `properties::NamedProperties`: arbitrary properties.
"""
mutable struct ProteinStructure{T} <: AbstractVector{ProteinChain{T}}
    name::String
    atoms::Vector{Atom{T}}
    chains::Vector{ProteinChain{T}}
    properties::NamedProperties

    function ProteinStructure{T}(name, atoms, chains, properties) where T
        return new{T}(name, atoms, chains, sortnames(namedproperties(properties)))
    end
end

ProteinStructure(name, atoms::Vector{Atom{T}}, chains, properties=(;)) where T =
    ProteinStructure{T}(name, atoms, chains, properties)

Base.convert(::Type{ProteinStructure{T}}, structure::ProteinStructure) where T =
    ProteinStructure(structure.name, convert(Vector{Atom{T}}, structure.atoms), convert(Vector{ProteinChain{T}}, structure.chains))

Base.size(structure::ProteinStructure) = (length(structure.chains),)

chainid_to_index(structure::ProteinStructure, id::AbstractString) = findfirst(c -> c.id == id, structure.chains)

Base.getindex(structure::ProteinStructure, i::Integer) = structure.chains[i]
Base.getindex(structure::ProteinStructure, id::AbstractString) = structure[chainid_to_index(structure, id)]

function Base.getindex(structure::ProteinStructure, inds::AbstractVector{<:Integer})
    properties = map(p -> p[inds], structure.properties)
    return ProteinStructure(structure.name, structure.atoms, map(i -> structure[i], inds), properties)
end

Base.setindex!(structure::ProteinStructure, chain::ProteinChain, i::Integer) = (structure.chains[i] = chain)
Base.setindex!(structure::ProteinStructure, chain::ProteinChain, id::AbstractString) = (structure.chains[chainid_to_index(structure, id)] = chain)

Base.summary(structure::ProteinStructure) = "$(length(structure))-chain $(typeof(structure)) \"$(structure.name)\""

function Base.show(io::IO, structure::ProteinStructure)
    print(io, "$(typeof(structure))(")
    for fieldname in fieldnames(ProteinStructure)
        show(io, getfield(structure, fieldname))
        fieldname != :properties && print(io, ", ")
    end
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", structure::ProteinStructure)
    print(io, summary(structure))
    for chain in structure
        print(io, "\n ", summary(chain))
    end
end

function map_chains!(f::Function, structure::ProteinStructure)
    for (i, chain) in enumerate(structure)
        structure[i] = f(chain)
    end
    return structure
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

Base.getproperty(structure::ProteinStructure, name::Symbol) =
    name in fieldnames(ProteinStructure) ? getfield(structure, name) : unpack(getfield(getfield(structure, :properties), name))

Base.propertynames(structure::ProteinStructure, private::Bool=false) = (setdiff(fieldnames(ProteinStructure), private ? () : (:properties,))..., propertynames(structure.properties)...)

function setproperties!(structure::ProteinStructure, properties::NamedTuple)
    structure.properties = setproperties(structure.properties, properties)
    structure
end

addproperties!(structure::ProteinStructure, properties::NamedTuple) =
    setproperties!(structure, addproperties(structure.properties, properties))

addproperties!(structure::ProteinStructure; properties...) =
    setproperties!(structure, addproperties(structure.properties, NamedTuple(properties)))

removeproperties!(structure::ProteinStructure, names::Symbol...) =
    setproperties!(structure, removeproperties(structure.properties, names...))
