"""
    ProteinStructure <: AbstractVector{ProteinChain}
"""
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

Base.summary(structure::ProteinStructure) = "$(length(structure))-chain ProteinStructure \"$(structure.name)\" with $(length(structure.properties)) properties"

function Base.show(io::IO, ::MIME"text/plain", structure::ProteinStructure)
    print(io, summary(structure), ":")
    n, i = 10, 0
    for chain in structure
        (i += 1) <= n || break
        print(io, "\n  ", summary(chain))
    end
    i < length(structure) && print("\n  …")
end
