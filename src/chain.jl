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

import Backboner

Backboner.Backbone(chain::ProteinChain) = Backbone(chain.backbone)
Backboner.ChainedBonds(chain::ProteinChain) = ChainedBonds(Backbone(chain))
Backboner.Frames(chain::ProteinChain, ideal_residue=STANDARD_RESIDUE) = Frames(Backbone(chain), ideal_residue)

# TODO: "contiguity" field/property to mask outputs with NaN where residues are not contiguous
psi_angles(chain::ProteinChain) = get_torsion_angles(Backbone(chain))[1:3:end]
omega_angles(chain::ProteinChain) = get_torsion_angles(Backbone(chain))[2:3:end]
phi_angles(chain::ProteinChain) = get_torsion_angles(Backbone(chain))[3:3:end]