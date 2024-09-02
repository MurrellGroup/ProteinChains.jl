@dynamic mutable struct ProteinStructure{T<:AbstractFloat} <: AbstractVector{ProteinChain{T}}
    name::String
    chains::Vector{ProteinChain{T}}
end

"""
    ProteinStructure{T<:AbstractFloat} <: AbstractVector{ProteinChain{T}}

## Examples

```jldoctest
julia> structure = pdb"1EYE"
1-chain ProteinStructure "1EYE.cif" with 2 dynamic properties:
  2 fields:
    name::String = "1EYE.cif"
    chains::Vector{ProteinChain{Float64}} = <field value exceeds max length>
  2 properties:
    ids::Vector{String} = ["A"]
    lengths::Vector{Int64} = [253]

julia> structure["A"] isa ProteinChain
true
```
"""
ProteinStructure

Base.size(structure::ProteinStructure) = (length(structure.chains),)
Base.axes(structure::ProteinStructure) = ([chain.id for chain in structure.chains],)
Base.getindex(structure::ProteinStructure, id::AbstractString) = structure.chains[findfirst(c -> c.id == id, structure.chains)]
Base.collect(structure::ProteinStructure) = structure.chains

Base.summary(structure::ProteinStructure) = "$(length(structure))-chain ProteinStructure \"$(structure.name)\""

function offset!(structure::ProteinStructure, coords::Vector{<:Real})
    @assert length(coords) == 3
    for chain in structure
        offset!(chain, coords)
    end
    return structure
end 