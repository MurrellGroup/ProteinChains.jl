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

function ProteinChain(
    chain::BioStructures.Chain;
    properties=(; ins_codes=BioStructures.inscode)
)
    residues = BioStructures.collectresidues(chain, backbone_residue_selector)
    return ProteinChain(
        BioStructures.chainid(chain),
        get_atoms(Atom{Float64}, residues),
        get_sequence(residues),
        map(BioStructures.resnumber, residues);
        (p => Indexable(map(f, residues)) for (p, f) in pairs(properties))...
    )
end

function ProteinStructure(model::BioStructures.Model)
    return ProteinStructure(
        model.structure.name,
        map(atom -> convert(Atom{Float64}, atom), BioStructures.collectatoms(BioStructures.collectresidues(model, !backbone_residue_selector))),
        [ProteinChain(chain) for chain in model]
    )
end

function ProteinStructure(struc::BioStructures.MolecularStructure)
    proteinstructure = ProteinStructure(first(BioStructures.collectmodels(struc)))
    return proteinstructure
end

function ProteinStructure(name::AbstractString, mmcifdict)
    return ProteinStructure(BioStructures.MolecularStructure(mmcifdict; structure_name=name))
end

function Base.read(
    path::AbstractString, ::Type{ProteinStructure}, ::Type{MMCIFFormat};
    mmcifdict=BioStructures.MMCIFDict(path),
)
    return ProteinStructure(basename(path), mmcifdict)
end

Base.read(path::AbstractString, P::Type{ProteinStructure}, ::Type{PDBFormat}; kwargs...) = P(read(path, PDBFormat); kwargs...)
Base.read(path::AbstractString, P::Type{ProteinStructure}; kwargs...) = read(path, P, get_format(path); kwargs...)

readcif(path::AbstractString) = read(path, ProteinStructure, MMCIFFormat)
readpdb(path::AbstractString) = read(path, ProteinStructure, PDBFormat)
