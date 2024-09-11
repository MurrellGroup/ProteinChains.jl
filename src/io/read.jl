backbone_atom_selector(atom::BioStructures.AbstractAtom) = BioStructures.atomnameselector(atom, BACKBONE_ATOM_NAMES)

backbone_residue_selector(residue::BioStructures.AbstractResidue) =
    oneletter_resname(residue) in AMINOACIDS &&
    BioStructures.countatoms(residue, backbone_atom_selector) == 3 &&
    BioStructures.standardselector(residue)

get_sequence(residues::Vector{BioStructures.AbstractResidue}) = join(map(oneletter_resname, residues))

get_atom(residue::BioStructures.AbstractResidue, name::AbstractString) =
    BioStructures.collectatoms(residue, a -> BioStructures.atomnameselector(a, [name])) |> only

get_backbone(residue::BioStructures.AbstractResidue) =
    stack(name -> BioStructures.coords(get_atom(residue, name)), BACKBONE_ATOM_NAMES)

get_backbone(residues::Vector{BioStructures.AbstractResidue}) =
    mapreduce(get_backbone, hcat, residues; init=Matrix{Float64}(undef, 3, 0))

Atom(atom::BioStructures.Atom) = Atom(atom.name, atom.element, atom.coords)
Atom(atom::BioStructures.DisorderedAtom) = Atom(BioStructures.defaultatom(atom))

function get_nonbackbone_atoms(residues::Vector{BioStructures.AbstractResidue})
    atoms = Vector{Atom{Float64}}[]
    for residue in residues
        residue = residue isa BioStructures.DisorderedResidue ? BioStructures.defaultresidue(residue) : residue
        residue_atoms = map(Atom, BioStructures.collectatoms(residue, !backbone_atom_selector))
        push!(atoms, residue_atoms)
    end
    return atoms
end

function ProteinChain(residues::Vector{BioStructures.AbstractResidue})
    id = only(unique(map(BioStructures.chainid, residues)))
    sequence = get_sequence(residues)
    backbone = reshape(get_backbone(residues), 3, 3, :)
    atoms = get_nonbackbone_atoms(residues)
    numbering = map(BioStructures.resnumber, residues)
    return ProteinChain(id, sequence, backbone, numbering, atoms)
end

function ProteinChain(chain::BioStructures.Chain, selector=backbone_residue_selector)
    residues = BioStructures.collectresidues(chain, selector)
    isempty(residues) && return ProteinChain(BioStructures.chainid(chain), "", zeros(3, 3, 0), Int[], Vector{Atom{Float64}}[])
    return ProteinChain(residues)
end

function ProteinStructure(struc::BioStructures.MolecularStructure, selector=backbone_residue_selector; path=nothing)
    chains = ProteinChain[]
    for model in struc, chain in model
        push!(chains, ProteinChain(chain, selector))
    end
    structure = ProteinStructure(struc.name, chains)
    !isnothing(path) && endswith(path, ".cif") && renumber!(structure, BioStructures.MMCIFDict(path))
    return structure
end

Base.read(path::AbstractString, ::Type{ProteinStructure}, format::Type{<:ProteinFileFormat}) = ProteinStructure(read(path, format); path)
Base.read(path::AbstractString, T::Type{ProteinStructure}) = read(path, T, get_format(path))

readcif(path::AbstractString) = read(path, ProteinStructure, MMCIFFormat)
readpdb(path::AbstractString) = read(path, ProteinStructure, PDBFormat)
