"""
    ProteinStructure <: AbstractVector{AbstractProteinChain}

```jldoctest
julia> pdb"1EYE"
[ Info: Downloading file from PDB: 1EYE
1-chain ProteinStructure "1EYE.cif"
 253-residue AnnotatedProteinChain{Float64} (A)
```
"""
mutable struct ProteinStructure <: AbstractVector{AbstractProteinChain}
    name::String
    chains::Vector{AbstractProteinChain}
end

Base.size(structure::ProteinStructure) = (length(structure.chains),)
Base.getindex(structure::ProteinStructure, i) = structure.chains[i]

Base.getindex(structure::ProteinStructure, id::AbstractString) = structure[findfirst(c -> c.id == id, structure.chains)]

Base.summary(structure::ProteinStructure) = "$(length(structure))-chain ProteinStructure \"$(structure.name)\""

function Base.show(io::IO, ::MIME"text/plain", structure::ProteinStructure)
    print(io, summary(structure))
    for chain in structure
        print(io, "\n ", summary(chain))
    end
end

function offset!(structure::ProteinStructure, coords::Vector{<:Real})
    for chain in structure
        offset!(chain, coords)
    end
    return structure
end 