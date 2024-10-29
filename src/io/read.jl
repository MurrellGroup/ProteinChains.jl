backbone_atom_selector(atom::BioStructures.AbstractAtom) = BioStructures.atomnameselector(atom, BACKBONE_ATOM_NAMES)

backbone_residue_selector(residue::BioStructures.AbstractResidue) =
    oneletter_resname(residue) in AMINOACIDS &&
    BioStructures.countatoms(residue, backbone_atom_selector) == 3 &&
    BioStructures.standardselector(residue)

get_sequence(residues::Vector{BioStructures.AbstractResidue}) = join(map(oneletter_resname, residues))

get_atom(residue::BioStructures.AbstractResidue, name::AbstractString) =
    BioStructures.collectatoms(residue, a -> BioStructures.atomnameselector(a, [name])) |> only

Base.convert(::Type{Atom{T}}, atom::BioStructures.Atom) where T = convert(Atom{T}, Atom(atom.name, atom.element, atom.coords))
Base.convert(::Type{Atom{T}}, atom::BioStructures.DisorderedAtom) where T = convert(Atom{T}, BioStructures.defaultatom(atom))

function get_atoms(::Type{Atom{T}}, residues::Vector{BioStructures.AbstractResidue}) where T
    atoms = Vector{Atom{T}}[]
    for residue in residues
        residue = residue isa BioStructures.DisorderedResidue ? BioStructures.defaultresidue(residue) : residue
        residue_atoms = map(atom -> convert(Atom{T}, atom), BioStructures.collectatoms(residue))
        push!(atoms, residue_atoms)
    end
    return atoms
end

function ProteinChain{T}(chain::BioStructures.Chain) where T
    residues = BioStructures.collectresidues(chain, backbone_residue_selector)
    proteinchain = if isempty(residues)
        ProteinChain(BioStructures.chainid(chain), Vector{Atom{T}}[], "", Int[])
    else
        id = BioStructures.chainid(chain)
        atoms = get_atoms(Atom{T}, residues)
        sequence = get_sequence(residues)
        numbering = map(BioStructures.resnumber, residues)
        ProteinChain(id, atoms, sequence, numbering)
    end
    return addproperties(proteinchain;
        ins_codes = IndexableProperty(map(Int8 ∘ BioStructures.inscode, residues)),
    )
end

function ProteinStructure{T}(model::BioStructures.Model) where T
    return ProteinStructure(
        model.structure.name,
        map(atom -> convert(Atom{T}, atom), BioStructures.collectatoms(BioStructures.collectresidues(model, !backbone_residue_selector))),
        [ProteinChain{T}(chain) for chain in model]
    )
end

function ProteinStructure{T}(struc::BioStructures.MolecularStructure) where T
    proteinstructure = ProteinStructure{T}(first(BioStructures.collectmodels(struc)))
    return proteinstructure
end

function ProteinStructure{T}(name::AbstractString, mmcifdict) where T
    return ProteinStructure{T}(BioStructures.MolecularStructure(mmcifdict; structure_name=name))
end

function Base.read(
    path::AbstractString, ::Type{ProteinStructure{T}}, ::Type{MMCIFFormat};
    mmcifdict=BioStructures.MMCIFDict(path),
) where T
    return ProteinStructure{T}(basename(path), mmcifdict)
end

Base.read(path::AbstractString, P::Type{ProteinStructure{T}}, ::Type{PDBFormat}; kwargs...) where T = P(read(path, PDBFormat); kwargs...)
Base.read(path::AbstractString, P::Type{<:ProteinStructure}; kwargs...) = read(path, P, get_format(path); kwargs...)
Base.read(path::AbstractString, ::Type{ProteinStructure}, args...; kwargs...) = read(path, ProteinStructure{Float64}, args...; kwargs...)

readcif(path::AbstractString) = read(path, ProteinStructure, MMCIFFormat)
readpdb(path::AbstractString) = read(path, ProteinStructure, PDBFormat)
