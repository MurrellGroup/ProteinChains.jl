"""
    AbstractProteinChain{T<:AbstractFloat}

An abstract type representing a protein chain, with a type parameter
`T` that specifies the floating-point type of the coordinates.

## Subtypes
- `ProteinChain{T}`: A concrete subtype of `AbstractProteinChain` that represents a protein chain with a spe.
"""
abstract type AbstractProteinChain{T<:AbstractFloat} end

Base.summary(chain::AbstractProteinChain) = "$(length(chain))-residue $(typeof(chain)) ($(chain.id))"
Base.length(chain::AbstractProteinChain) = size(chain.backbone, 3)

mutable struct ProteinChain{T} <: AbstractProteinChain{T}
    id::String
    sequence::String
    backbone::Array{T,3}
    numbering::Vector{Int}
    atoms::Vector{Vector{Atom{T}}}

    function ProteinChain(
        id::AbstractString, sequence::AbstractString, backbone::AbstractArray{T,3},
        numbering::AbstractVector{<:Integer} = collect(1:length(sequence)),
        atoms::Vector{Vector{Atom{T}}} = [Atom{T}[] for _ in sequence],
    ) where T
        @assert length(sequence) == size(backbone, 3) == length(numbering) == length(atoms)
        new{T}(id, sequence, backbone, numbering, atoms)
    end
end

"""
    ProteinChain{T<:AbstractFloat}
"""
ProteinChain

Base.getindex(chain::ProteinChain, i::AbstractVector{<:Integer}) =
    ProteinChain(chain.id, chain.sequence[i], chain.backbone[:,:,i], chain.numbering[i], chain.atoms[i])

import Backboner

Backboner.Backbone(chain::AbstractProteinChain) = Backbone(chain.backbone)
Backboner.ChainedBonds(chain::AbstractProteinChain) = ChainedBonds(Backbone(chain))
Backboner.Frames(chain::AbstractProteinChain, ideal_residue=STANDARD_RESIDUE) = Frames(Backbone(chain), ideal_residue)

psi_angles(chain::AbstractProteinChain) = get_torsion_angles(Backbone(chain))[1:3:end]
omega_angles(chain::AbstractProteinChain) = get_torsion_angles(Backbone(chain))[2:3:end]
phi_angles(chain::AbstractProteinChain) = get_torsion_angles(Backbone(chain))[3:3:end]

function offset!(chain::AbstractProteinChain, coords::Vector{<:Real})
    chain.backbone .+= coords
    for residue_atoms in chain.atoms
        for atom in residue_atoms
            offset!(atom, coords)
        end
    end
    return chain
end
