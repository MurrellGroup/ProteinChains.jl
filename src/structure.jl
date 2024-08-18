"""
    ProteinStructure{T<:AbstractFloat} <: AbstractVector{ProteinChain{T}}
"""
mutable struct ProteinStructure{T<:AbstractFloat} <: AbstractVector{ProteinChain{T}}
    name::String
    chains::Vector{ProteinChain{T}}
    properties::Properties
end

function ProteinStructure(name::String, chains::Vector{ProteinChain{T}}; kwargs...) where T
    return ProteinStructure{T}(name, chains, Properties(;
        ids=map(chain -> chain.id, chains),
        lengths=map(countresidues, chains),
        kwargs...))
end

Base.hasproperty(structure::ProteinStructure, name::Symbol) = hasproperty(HasProperties(structure), name)
Base.getproperty(structure::ProteinStructure, name::Symbol) = getproperty(HasProperties(structure), name)
Base.setproperty!(structure::ProteinStructure, name::Symbol, value) = setproperty!(HasProperties(structure), name, value)
Base.propertynames(structure::ProteinStructure) = propertynames(HasProperties(structure))

Base.size(structure::ProteinStructure) = (length(structure.chains),)
Base.getindex(structure::ProteinStructure, i) = structure.chains[i]

Base.getindex(structure::ProteinStructure, id::AbstractString) = structure[findfirst(c -> c.id == id, structure.chains)]

Base.summary(structure::ProteinStructure) = "$(length(structure))-chain ProteinStructure \"$(structure.name)\" with $(length(structure.properties)) properties"

Base.show(io::IO, ::MIME"text/plain", structure::ProteinStructure) = showproperties(io, structure)