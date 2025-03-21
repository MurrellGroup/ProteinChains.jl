"""
    ProteinChain{T<:Real}

Represents a protein chain with a basic set of fields from which some other properties might be derived.
"""
@dynamic mutable struct ProteinChain{T<:Real}
    id::String
    atoms::Vector{Vector{Atom{T}}}
    sequence::String
    numbering::Vector{Int}
end

ProteinChain(id, atoms, sequence; kws...) = ProteinChain(id, atoms, sequence, collect(1:length(atoms)); kws...)

function Base.convert(::Type{ProteinChain{T}}, chain::ProteinChain) where T
    ProteinChain(
        chain.id,
        convert(Vector{Vector{Atom{T}}}, chain.atoms),
        chain.sequence,
        chain.numbering;
        propertypairs(chain, NoFields)...)
end

Base.length(chain::ProteinChain) = length(chain.atoms)

function Base.getindex(chain::ProteinChain, i)
    args = chain.id, chain.atoms[i], chain.sequence[i], chain.numbering[i]
    kws = Iterators.map(propertypairs(chain, NoFields)) do (name, value)
        name => if value isa AbstractProperty
            checkproperty(chain, value)
            value isa Indexable ? value[i] : value
        else
            value
        end
    end
    return ProteinChain(args...; kws...)
end

Base.summary(chain::ProteinChain) = "$(length(chain))-residue $(typeof(chain)) ($(chain.id))"

function Base.show(io::IO, ::MIME"text/plain", chain::ProteinChain)
    print(io, summary(chain))
end

function get_atoms(backbone_coords::Array{T,3}) where T
    @assert size(backbone_coords)[1:2] == (3,3)
    atoms = [Vector{Atom{T}}(undef, 3) for _ in 1:size(backbone_coords, 3)]
    for (i, slice) in enumerate(eachslice(backbone_coords, dims=3))
        for (j, pos) in enumerate(eachcol(slice))
            atoms[i][j] = Atom(BACKBONE_ATOM_NAMES[j], BACKBONE_ATOM_SYMBOLS[j], pos)
        end
    end
    return atoms
end

get_atoms(backbone::Backbone) = get_atoms(reshape(backbone.coords, 3, 3, :))
get_atoms(chain::ProteinChain) = chain.atoms

function get_backbone(atoms::Vector{Vector{Atom{T}}}) where T
    backbone_coords = Array{T,3}(undef, 3, 3, length(atoms))
    encoded_names = encode_atom_name.(BACKBONE_ATOM_NAMES, BACKBONE_ATOM_SYMBOLS)
    for (i, residue_atoms) in enumerate(atoms)
        for (j, name) in enumerate(encoded_names)
            backbone_coords[:,j,i] = atom_coords(argmax(atom -> atom.name == name, residue_atoms))
        end
    end
    return backbone_coords
end

get_backbone(chain::ProteinChain) = get_backbone(chain.atoms)

Backboner.Backbone(chain::ProteinChain) = Backbone(get_backbone(chain))
Backboner.ChainedBonds(chain::ProteinChain) = ChainedBonds(Backbone(chain))
Backboner.Frames(chain::ProteinChain, ideal_residue=STANDARD_RESIDUE) = Frames(Backbone(chain), ideal_residue)

Backboner.get_bond_angles(chain::ProteinChain) = get_bond_angles(Backbone(chain))
Backboner.get_bond_lengths(chain::ProteinChain) = get_bond_lengths(Backbone(chain))
Backboner.get_torsion_angles(chain::ProteinChain) = get_torsion_angles(Backbone(chain))

psi_angles(chain::ProteinChain) = get_torsion_angles(Backbone(chain))[1:3:end]
omega_angles(chain::ProteinChain) = get_torsion_angles(Backbone(chain))[2:3:end]
phi_angles(chain::ProteinChain) = get_torsion_angles(Backbone(chain))[3:3:end]

function map_atoms!(f::Function, chain::ProteinChain, args...)
    for residue_atoms in chain.atoms
        for i in eachindex(residue_atoms)
            residue_atoms[i] = f(residue_atoms[i], args...)
        end
    end
    return chain
end
