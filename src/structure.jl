@properties mutable struct ProteinStructure{T<:AbstractFloat} <: AbstractVector{ProteinChain{T}}
    name::String
    chains::Vector{ProteinChain{T}}
end

"""
    ProteinStructure{T<:AbstractFloat} <: AbstractVector{ProteinChain{T}}
"""
ProteinStructure

function ProteinStructure(name::String, chains::Vector{ProteinChain{T}}; kwargs...) where T
    return ProteinStructure{T}(name, chains, Properties(;
        ids=map(chain -> chain.id, chains),
        lengths=map(countresidues, chains),
        kwargs...))
end

Base.size(structure::ProteinStructure) = (length(structure.chains),)
Base.getindex(structure::ProteinStructure, i) = structure.chains[i]

Base.getindex(structure::ProteinStructure, id::AbstractString) = structure[findfirst(c -> c.id == id, structure.chains)]

Base.summary(structure::ProteinStructure) = "$(length(structure))-chain ProteinStructure \"$(structure.name)\" with $(length(structure.properties)) properties"

Base.show(io::IO, ::MIME"text/plain", structure::ProteinStructure) = showproperties(io, structure)

function offset!(structure::ProteinStructure, coords::Vector{<:Real})
    @assert length(coords) == 3
    for chain in structure
        offset!(chain, coords)
    end
    return structure
end 