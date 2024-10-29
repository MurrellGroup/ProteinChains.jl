"""
    ProteinChain{T<:Real}

Represents a protein chain with a basic set of fields from which some other properties might be derived.
The [`addproperties`](@ref) function can be used to instantiate new chains with additional properties.

## Fields
- `id::String`: Identifier for the protein chain.
- `atoms::Vector{Vector{Atom{T}}}`: List of atoms in each residue.
- `sequence::String`: Amino acid sequence of the protein.
- `numbering::Vector{Int32}`: Residue numbering (author). See [`renumber`](@ref) for renumbering.
- `properties::Ps`: Named properties associated with the chain.

See also [`addproperties`](@ref), [`StandardProperty`](@ref), [`IndexableProperty`](@ref).
```
"""
struct ProteinChain{T<:Real}
    id::String
    atoms::Vector{Vector{Atom{T}}}
    sequence::String
    numbering::Vector{Int32}
    properties::NamedProperties

    function ProteinChain{T}(
        id::AbstractString,
        atoms::Vector{Vector{Atom{T}}},
        sequence::AbstractString,
        numbering::Vector{Int32},
        properties::NamedProperties,
    ) where T
        @assert length(atoms) == sizeof(sequence) == length(numbering)
        chain = new{T}(id, atoms, sequence, numbering, sortnames(properties))
        for property in properties
            checkproperty(chain, property)
        end
        return chain
    end
end

function ProteinChain(id, atoms::Vector{Vector{Atom{T}}}, sequence, numbering::Vector{<:Integer}, properties::NamedTuple) where T
    ProteinChain{T}(id, atoms, sequence, Int32.(numbering), namedproperties(properties))
end

ProteinChain(id, atoms, sequence, numbering) = ProteinChain(id, atoms, sequence, numbering, (;))

Base.convert(::Type{ProteinChain{T}}, chain::ProteinChain) where T =
    ProteinChain(chain.id, convert(Vector{Vector{Atom{T}}}, chain.atoms), chain.sequence, chain.numbering, chain.properties)

function Base.:(==)(chain1::ProteinChain, chain2::ProteinChain)
    propertynames(chain1, false) != propertynames(chain2, false) && return false
    !any(getproperty(chain1, name) != getproperty(chain2, name) for name in propertynames(chain1, false))
end

Base.length(chain::ProteinChain) = length(chain.atoms)

Base.getproperty(chain::ProteinChain, name::Symbol) =
    name in fieldnames(ProteinChain) ? getfield(chain, name) : unpack(getfield(getfield(chain, :properties), name))

Base.propertynames(chain::ProteinChain, private::Bool=false) = (setdiff(fieldnames(ProteinChain), private ? () : (:properties,))..., propertynames(chain.properties)...)

function Base.getindex(chain::ProteinChain, i::Union{AbstractVector,Colon})
    properties = map(p -> p[i], chain.properties)
    ProteinChain(chain.id, chain.atoms[i], chain.sequence[i], chain.numbering[i], properties)
end

setproperties(chain::ProteinChain, ps::NamedTuple) = ProteinChain(chain.id, chain.atoms, chain.sequence, chain.numbering, ps)

"""
    addproperties(chain::ProteinChain; properties...)

Creates a new `ProteinChain` instance with the added properties.
Indexing of property values can be specified with a wrapper type,
such as `IndexableProperty`.    

See also [`removeproperties`](@ref), [`IndexableProperty`](@ref).
"""
addproperties(chain::ProteinChain; properties...) = setproperties(chain, merge(chain.properties, properties))

"""
    removeproperties(chain::ProteinChain, names::Symbol...)

Creates a new `ProteinChain` instance with the property names in `names` removed.

See also [`addproperties`](@ref)
"""
removeproperties(chain::ProteinChain, names::Symbol...) = setproperties(chain, NamedTuple{filter(name -> name ∉ names, propertynames(chain.properties))}(chain.properties))

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

"""
    addproperties(chain::ProteinChain, names::Symbol...)
    addproperties(chain::ProteinStructure, names::Symbol...)

Add predefined properties to a chain or chains of a structure.
"""
addproperties(chain::ProteinChain, names::Symbol...) = addproperties(chain; NamedTuple{names}(calculate_property(chain, name) for name in names)...)

calculate_property(x, name::Symbol, args...) = calculate_property(x, Val(name), args...)

calculate_property(chain::ProteinChain, ::Val{:ideal_residue}) = collect(STANDARD_RESIDUE)
calculate_property(chain::ProteinChain, ::Val{:bond_lengths}) = Backboner.get_bond_lengths(Backboner.Backbone(get_backbone(chain)))
calculate_property(chain::ProteinChain, ::Val{:bond_angles}) = Backboner.get_bond_angles(Backboner.Backbone(get_backbone(chain)))
calculate_property(chain::ProteinChain, ::Val{:torsion_angles}) = Backboner.get_torsion_angles(Backboner.Backbone(get_backbone(chain)))
calculate_property(chain::ProteinChain, ::Val{:is_knotted}) = Backboner.is_knotted(Backboner.Backbone(get_backbone(chain)[:,2,:]))

calculate_property(chain::ProteinChain, ::Val{:backbone}) = get_backbone(chain) |> IndexableProperty
calculate_property(chain::ProteinChain, ::Val{:secondary_structure}) = Int8.(assign_secondary_structure(get_backbone(chain))) |> IndexableProperty
calculate_property(chain::ProteinChain, ::Val{:residue_rotations}) = Backboner.Frames(Backbone(get_backbone(chain)), hasproperty(chain, :ideal_residue) ? chain.ideal_residue : STANDARD_RESIDUE).rotations |> IndexableProperty
calculate_property(chain::ProteinChain, ::Val{:residue_translations}) = dropdims(Backboner.centroid(get_backbone(chain); dims=2); dims=2) |> IndexableProperty
calculate_property(chain::ProteinChain{T}, ::Val{:residue_torsion_angles}) where T = [reshape(calculate_property(chain, :torsion_angles)[], 3, :) fill(T(NaN), 3)] |> IndexableProperty
