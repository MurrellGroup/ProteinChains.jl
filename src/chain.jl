@dynamic mutable struct ProteinChain{T<:AbstractFloat}
    id::String
    sequence::String
    backbone::Array{T,3}
    atoms::Vector{Vector{Atom{T}}}
end

"""
    ProteinChain{T<:AbstractFloat}

## Examples

```jldoctest
julia> structure = pdb"1EYE";
[ Info: Downloading file from PDB: 1EYE

julia> structure[1]
253-residue ProteinChain "A":
  4 fields:
    id::String = "A"
    sequence::String = <field value exceeds max length>
    backbone::Array{Float64,3} = <field value exceeds max length>
    atoms::Vector{Vector{ProteinChains.Atom{Float64}}} = <field value exceeds max length>
  2 properties:
    numbering::Vector{Int64} = [5, 6, 7, 8, 9, 10, 11, 12, 13, 14  …  265, 266, 267, 268, 269, 270, 271, 272, 273, 274]
    modelnum::Int64 = 1
```
"""
ProteinChain

countresidues(chain::ProteinChain) = length(chain.sequence)

Base.summary(chain::ProteinChain) = "$(countresidues(chain))-residue ProteinChain \"$(chain.id)\""

function offset!(chain::ProteinChain, coords::Vector{<:Real})
    @assert length(coords) == 3
    chain.backbone .+= coords
    for residue_atoms in chain.atoms
        for atom in residue_atoms
            offset!(atom, coords)
        end
    end
    return chain
end 