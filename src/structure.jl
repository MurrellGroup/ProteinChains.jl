mutable struct ProteinStructure <: AbstractVector{ProteinChain}
    name::String
    chains::Vector{ProteinChain}
    properties::Dict{Symbol,Any}
end

function ProteinStructure(name::String, chains::Vector{ProteinChain}; kwargs...)
    ps = ProteinStructure(name, chains, Dict{Symbol,Any}())
    !isempty(kwargs) && push!(pc.properties, kwargs...)
    return ps
end

Base.getproperty(structure::ProteinStructure, property::Symbol) = _getproperty(structure, property)
Base.setproperty!(structure::ProteinStructure, property::Symbol, value) = _setproperty!(structure, property, value)

Base.size(structure::ProteinStructure) = (length(structure.chains),)
Base.getindex(structure::ProteinStructure, i) = structure.chains[i]

Base.getindex(structure::ProteinStructure, id::AbstractString) = structure[findfirst(c -> c.id == id, structure.chains)]