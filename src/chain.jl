@properties mutable struct ProteinChain{T<:AbstractFloat}
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
253-residue ProteinChain "A" with 2 properties:
  fields:
    id::String = "A"
    sequence::String = "PVQVMGVLNVTDDSFSDGGCYLDLDDAVKHGLAMAAAGAGIVDVGGETSRVIPVVKELAAQGITVIDTMRADVARAALQNGAQMVNDVSGGRADPAMG…
    backbone::Array{Float64,3} = [45.592 44.171 43.719; -10.864 -10.936 -9.688; 30.192 30.504 31.278;;; 42.568 42.02 40.707; -9.163 …
    atoms::Vector{Vector{ProteinChains.Atom{Float64}}} = Vector{Atom{Float64}}[[Atom{Float64}(0x4f202020, 0x08, 44.41, -9.176, 32.157), Atom{Float64}(0x4243…
  properties:
    numbering::Vector{Int64} = [5, 6, 7, 8, 9, 10, 11, 12, 13, 14  …  265, 266, 267, 268, 269, 270, 271, 272, 273, 274]
    modelnum::Int64 = 1
```
"""
ProteinChain

countresidues(chain::ProteinChain) = length(chain.sequence)

function ProteinChain(id::String, sequence::String, backbone::Array{T,3}, atoms::Vector{Vector{Atom{T}}}, properties::Properties) where T
    @assert size(backbone)[1:2] == (3, 3)
    return ProteinChain{T}(id, sequence, backbone, atoms, properties)
end

function ProteinChain(id::String, sequence::String, backbone::Array{T,3}, atoms::Vector{Vector{Atom{T}}}=Vector{Atom{T}}[]; kwargs...) where T
    return ProteinChain(id, sequence, backbone, atoms, Properties(; kwargs...))
end


Base.summary(chain::ProteinChain) = "$(countresidues(chain))-residue ProteinChain \"$(chain.id)\" with $(length(chain.properties)) properties"

Base.show(io::IO, ::MIME"text/plain", chain::ProteinChain) = showproperties(io, chain)

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