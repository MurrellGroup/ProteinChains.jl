"""
    ProteinChain{T<:Real,Ps<:NamedProperties}
"""
struct ProteinChain{T<:Real,Ps<:NamedProperties}
    id::String
    atoms::Vector{Vector{Atom{T}}}
    sequence::String
    numbering::Vector{Int32}
    properties::Ps
end

function ProteinChain(
    id, atoms::Vector{Vector{Atom{T}}}, sequence::String, numbering::Vector{<:Integer}, properties::Ps,
) where {T,Ps<:NamedProperties}
    len = length(atoms)
    @assert sizeof(sequence) == len
    @assert length(numbering) == len
    for property in properties
        property isa ResidueProperty && @assert size(property[], ndims(property[])) == len
    end
    ProteinChain{T,Ps}(id, atoms, sequence, Int32.(numbering), properties)
end

ProteinChain(id, atoms, sequence, numbering) = ProteinChain(id, atoms, sequence, numbering, (;))

Base.convert(::Type{ProteinChain{T}}, chain::ProteinChain) where T =
    ProteinChain(chain.id, convert(Vector{Vector{Atom{T}}}, chain.atoms), chain.sequence, chain.numbering, chain.properties)

Base.:(==)(chain1::ProteinChain, chain2::ProteinChain) = propertynames(chain1) == propertynames(chain2) &&
    !any(getproperty(chain1, name) != getproperty(chain2, name) for name in propertynames(chain1, false))

Base.length(chain::ProteinChain) = length(chain.atoms)

Base.getproperty(chain::ProteinChain, name::Symbol) =
    name in fieldnames(ProteinChain) ? getfield(chain, name) : getfield(getfield(chain, :properties), name)[]

Base.propertynames(chain::ProteinChain, private::Bool=false) = (setdiff(fieldnames(ProteinChain), private ? () : (:properties,))..., propertynames(chain.properties)...)

function Base.getindex(chain::ProteinChain, i::AbstractVector)
    properties = map(p -> p[i], chain.properties)
    ProteinChain(chain.id, chain.atoms[i], chain.sequence[i], chain.numbering[i], properties)
end

annotate(chain::ProteinChain; annotations...) =
    ProteinChain(chain.id, chain.atoms, chain.sequence, chain.numbering, merge(chain.properties, annotations))

Base.summary(chain::ProteinChain) = "$(length(chain))-residue $(typeof(chain)) ($(chain.id))"

# wrap io with IOContext(io, :compact=>false) to make parseable
function Base.show(io::IO, chain::ProteinChain)
    print(io, "ProteinChain(")
    for fieldname in fieldnames(ProteinChain)
        show(io, getproperty(chain, fieldname))
        fieldname == :properties || print(io, ", ")
    end
    print(io, ")")
end

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

get_backbone(chain::ProteinChain) = hasproperty(chain, :backbone) ? chain.backbone : get_backbone(chain.atoms)

Backboner.Backbone(chain::ProteinChain) = Backbone(get_backbone(chain))
Backboner.ChainedBonds(chain::ProteinChain) = ChainedBonds(Backbone(chain))
Backboner.Frames(chain::ProteinChain, ideal_residue=STANDARD_RESIDUE) = Frames(Backbone(chain), ideal_residue)

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

using AssigningSecondaryStructure: assign_secondary_structure

calculate_property(x, name::Symbol, args...) = calculate_property(x, Val(name), args...)
annotate(chain::ProteinChain, names::Vararg{Symbol}) = annotate(chain; NamedTuple{names}(calculate_property(chain, name) for name in names)...)

calculate_property(chain::ProteinChain, ::Val{:ideal_residue}) = collect(STANDARD_RESIDUE) |> ChainProperty
calculate_property(chain::ProteinChain, ::Val{:bond_lengths}) = Backboner.get_bond_lengths(Backboner.Backbone(get_backbone(chain))) |> ChainProperty
calculate_property(chain::ProteinChain, ::Val{:bond_angles}) = Backboner.get_bond_angles(Backboner.Backbone(get_backbone(chain))) |> ChainProperty
calculate_property(chain::ProteinChain, ::Val{:torsion_angles}) = Backboner.get_torsion_angles(Backboner.Backbone(get_backbone(chain))) |> ChainProperty
calculate_property(chain::ProteinChain, ::Val{:is_knotted}) = Backboner.is_knotted(Backboner.Backbone(get_backbone(chain)[:,2,:])) |> ChainProperty

calculate_property(chain::ProteinChain, ::Val{:backbone}) = get_backbone(chain) |> ResidueProperty
calculate_property(chain::ProteinChain, ::Val{:secondary_structure}) = Int8.(assign_secondary_structure(get_backbone(chain))) |> ResidueProperty
calculate_property(chain::ProteinChain, ::Val{:residue_rotations}) = Backboner.Frames(Backbone(get_backbone(chain)), hasproperty(chain, :ideal_residue) ? chain.ideal_residue : STANDARD_RESIDUE).rotations |> ResidueProperty
calculate_property(chain::ProteinChain, ::Val{:residue_translations}) = dropdims(Backboner.centroid(get_backbone(chain); dims=2); dims=2) |> ResidueProperty
calculate_property(chain::ProteinChain{T}, ::Val{:residue_torsion_angles}) where T = [reshape(calculate_property(chain, :torsion_angles)[], 3, :) fill(T(NaN), 3)] |> ResidueProperty
