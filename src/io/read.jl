backbone_atom_selector(atom::BioStructures.AbstractAtom) = BioStructures.atomnameselector(atom, BACKBONE_ATOM_NAMES)

backbone_residue_selector(residue::BioStructures.AbstractResidue) =
    oneletter_resname(residue) in AMINOACIDS &&
    BioStructures.countatoms(residue, backbone_atom_selector) == 3 &&
    BioStructures.standardselector(residue)

get_sequence(residues::Vector{BioStructures.AbstractResidue}) = join(map(oneletter_resname, residues))

get_atom(residue::BioStructures.AbstractResidue, name::AbstractString) =
    BioStructures.collectatoms(residue, a -> BioStructures.atomnameselector(a, [name])) |> only

Atom(atom::BioStructures.Atom) = Atom(atom.name, atom.element, atom.coords)
Atom(atom::BioStructures.DisorderedAtom) = Atom(BioStructures.defaultatom(atom))

function get_atoms(residues::Vector{BioStructures.AbstractResidue})
    atoms = Vector{Atom{Float64}}[]
    for residue in residues
        residue = residue isa BioStructures.DisorderedResidue ? BioStructures.defaultresidue(residue) : residue
        residue_atoms = map(Atom, BioStructures.collectatoms(residue))
        push!(atoms, residue_atoms)
    end
    return atoms
end

function ProteinChain(residues::Vector{BioStructures.AbstractResidue})
    id = only(unique(map(BioStructures.chainid, residues)))
    atoms = get_atoms(residues)
    sequence = get_sequence(residues)
    numbering = map(BioStructures.resnumber, residues)
    return ProteinChain(id, atoms, sequence, numbering)
end

function ProteinChain(chain::BioStructures.Chain, selector=backbone_residue_selector)
    residues = BioStructures.collectresidues(chain, selector)
    isempty(residues) && return ProteinChain(BioStructures.chainid(chain), "", zeros(3, 3, 0), Int[], Vector{Atom{Float64}}[])
    return ProteinChain(residues)
end

function ProteinStructure(struc::BioStructures.MolecularStructure, selector=backbone_residue_selector; mmcifdict=nothing)
    proteinchains = ProteinChain{Float64}[]
    atoms = Atom.(BioStructures.collectatoms(BioStructures.collectresidues(struc, !selector)))
    for model in struc, chain in model
        push!(proteinchains, ProteinChain(chain, selector))
    end
    proteinstructure = ProteinStructure(struc.name, atoms, proteinchains)
    mmcifdict isa BioStructures.MMCIFDict && renumber!(proteinstructure, mmcifdict)
    return proteinstructure
end

Base.read(path::AbstractString, ::Type{ProteinStructure}, format::Type{MMCIFFormat}) = ProteinStructure(read(path, format); mmcifdict=BioStructures.MMCIFDict(path))
Base.read(path::AbstractString, ::Type{ProteinStructure}, format::Type{<:ProteinFileFormat}) = ProteinStructure(read(path, format))
Base.read(path::AbstractString, T::Type{ProteinStructure}) = read(path, T, get_format(path))

readcif(path::AbstractString) = read(path, ProteinStructure, MMCIFFormat)
readpdb(path::AbstractString) = read(path, ProteinStructure, PDBFormat)
