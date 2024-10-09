"""
    ProteinStructure{T} <: AbstractVector{ProteinChain{T}}
"""
struct ProteinStructure{T} <: AbstractVector{ProteinChain{T}}
    name::String
    atoms::Vector{Atom{T}}
    chains::Vector{<:ProteinChain{T}}
    numbering::Vector{Int}
end

ProteinStructure(name, atoms, chains) = ProteinStructure(name, atoms, chains, collect(1:length(chains)))

Base.convert(::Type{ProteinStructure{T}}, structure::ProteinStructure) where T =
    ProteinStructure(structure.name, convert(Vector{Atom{T}}, structure.atoms), convert(Vector{ProteinChain{T}}, structure.chains), structure.numbering)

Base.size(structure::ProteinStructure) = (length(structure.chains),)

Base.getindex(structure::ProteinStructure, i::Integer) = structure.chains[i]

function Base.getindex(structure::ProteinStructure, i::AbstractVector)
    ProteinStructure(structure.name, structure.atoms, structure.chains[i], structure.numbering[i])
end

Base.getindex(structure::ProteinStructure, id::AbstractString) = structure[findfirst(c -> c.id == id, structure.chains)]

Base.summary(structure::ProteinStructure) = "$(length(structure))-chain $(typeof(structure)) \"$(structure.name)\""

function Base.show(io::IO, structure::ProteinStructure)
    print(io, "$(typeof(structure))(")
    for fieldname in fieldnames(ProteinStructure)
        show(io, getfield(structure, fieldname))
        fieldname != :numbering && print(io, ", ")
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

annotate(structure::ProteinStructure, names::Vararg{Symbol}) =
    ProteinStructure(structure.name, structure.atoms, annotate.(structure.chains, names...), structure.numbering)