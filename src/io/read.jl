backbone_atom_selector(atom::BioStructures.AbstractAtom) = BioStructures.atomnameselector(atom, BACKBONE_ATOM_NAMES)
backbone_residue_selector(residue::BioStructures.AbstractResidue) =
    oneletter_resname(residue) in AMINOACIDS &&
    BioStructures.countatoms(residue, backbone_atom_selector) == 3 &&
    BioStructures.standardselector(residue) &&
    !BioStructures.disorderselector(residue) &&
    !any(atom -> atom isa BioStructures.DisorderedAtom, BioStructures.collectatoms(residue))

function get_atom(residue::BioStructures.AbstractResidue, atom_name::AbstractString)
    selector = atom -> BioStructures.atomnameselector(atom, [atom_name])
    residue_atoms = BioStructures.collectatoms(residue, selector)
    return only(residue_atoms)
end

function get_backbone(residue::BioStructures.AbstractResidue)
    atom_name_to_atom_coords = atom_name -> BioStructures.coords(get_atom(residue, atom_name))
    residue_coords = stack(atom_name_to_atom_coords, BACKBONE_ATOM_NAMES)
    return Backbone(residue_coords)
end

function get_backbone(residues::Vector{BioStructures.AbstractResidue})
    chain_coords = mapreduce(residue -> get_backbone(residue).coords, hcat, residues; init=Matrix{Float64}(undef, 3, 0))
    return Backbone(chain_coords)
end

get_backbone(chain::BioStructures.Chain, selector) = get_backbone(BioStructures.collectresidues(chain, selector))

aminoacid_sequence(residues::Vector{BioStructures.AbstractResidue}) = join(oneletter_resname.(residues))
aminoacid_sequence(chain::BioStructures.Chain, selector) = aminoacid_sequence(BioStructures.collectresidues(chain, selector))

function get_nonbackbone_atoms(residues::Vector{BioStructures.AbstractResidue})
    atoms = Vector{Atom{Float64}}[]
    for residue in residues
        residue_atoms_bs = BioStructures.collectatoms(residue, a -> !backbone_atom_selector(a) && BioStructures.standardselector(a) && !BioStructures.disorderselector(a))
        @assert !isempty(residue_atoms_bs)
        push!(atoms, map(atom_bs -> Atom(atom_bs.name, atom_bs.element, atom_bs.coords), residue_atoms_bs))
    end
    return atoms
end

function ProteinChain(residues::Vector{BioStructures.AbstractResidue})
    id = only(unique(BioStructures.chainid.(residues)))
    sequence = aminoacid_sequence(residues)
    backbone = reshape(get_backbone(residues).coords, 3, 3, :)
    atoms = get_nonbackbone_atoms(residues)
    numbering = BioStructures.resnumber.(residues)
    modelnum = only(unique(BioStructures.modelnumber.(residues)))
    #ins_codes = BioStructures.inscode.(residues)
    return ProteinChain(id, sequence, backbone, atoms; numbering, modelnum)
end

function ProteinChain(chain::BioStructures.Chain, selector=backbone_residue_selector)
    residues = BioStructures.collectresidues(chain, selector)
    isempty(residues) && return ProteinChain(BioStructures.chainid(chain), "", zeros(3, 3, 0), Int[])
    return ProteinChain(residues)
end

function ProteinStructure(struc::BioStructures.MolecularStructure, selector=backbone_residue_selector)
    name = struc.name
    chains = ProteinChain{Float64}[]
    for model in struc, chain in model
        if !isempty(chain)
            pc = ProteinChain(chain, selector)
            countresidues(pc) > 0 && push!(chains, pc)
        end
    end
    ProteinStructure(name, chains)
end
