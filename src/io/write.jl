function BioStructures.Chain(proteinchain::ProteinChain, model::BioStructures.Model)
    numbering = hasproperty(proteinchain, :numbering) ? proteinchain.numbering : collect(1:countresidues(proteinchain))
    residue_list = Vector{String}()
    residues = Dict{String, BioStructures.AbstractResidue}()
    chain = BioStructures.Chain(proteinchain.id, residue_list, residues, model)
    atom_serial = 0
    for residue_index in 1:countresidues(proteinchain)
        atom_list = Vector{String}()
        atoms = Dict{String, BioStructures.AbstractAtom}()
        resname = threeletter_resname(proteinchain.sequence[residue_index])
        number = numbering[residue_index]
        residue = BioStructures.Residue(resname, number, ' ', false, atom_list, atoms, chain, '-') # TODO: secondary structure
        for (i, (atom_name, element)) in enumerate(zip(BACKBONE_ATOM_NAMES, BACKBONE_ATOM_SYMBOLS))
            atom_serial += 1
            coords = proteinchain.backbone[:, i, residue_index]
            atom = BioStructures.Atom(atom_serial, atom_name, ' ', coords, 1.0, 0.0, element, "  ", residue)
            push!(atom_list, atom.name)
            atoms[atom.name] = atom
        end
        for atom in proteinchain.atoms[residue_index]
            atom_serial += 1
            atom_name = decode_atom_name(atom.atom_name)
            coords = [atom.x, atom.y, atom.z]
            element = atomic_number_to_element_symbol(atom.atomic_number)
            atom = BioStructures.Atom(atom_serial, atom_name, ' ', coords, 1.0, 0.0, element, "  ", residue)
            push!(atom_list, atom.name)
            atoms[atom.name] = atom
        end
        key = string(residue.number)
        push!(residue_list, key)
        residues[key] = residue
    end
    return chain
end

function BioStructures.MolecularStructure(proteinstruc::ProteinStructure)
    models = Dict{Int64, BioStructures.Model}()
    struc = BioStructures.MolecularStructure(proteinstruc.name, models)
    for proteinchain in proteinstruc
        modelnum = hasproperty(proteinchain, :modelnum) ? proteinchain.modelnum : 1
        modelnum in keys(models) || (models[modelnum] = BioStructures.Model(modelnum, Dict{String, BioStructures.Chain}(), struc))
        model = models[modelnum]
        model.chains[proteinchain.id] = BioStructures.Chain(proteinchain, model)
    end
    return struc
end

function writechains(path::AbstractString, proteinstruc::ProteinStructure, format::Type{<:ProteinFileFormat})
    struc = BioStructures.MolecularStructure(proteinstruc)
    if format == MMCIFFormat
        BioStructures.writemmcif(path, struc)
    elseif format == PDBFormat
        BioStructures.writepdb(path, struc)
    else
        error("Unsupported format: $format")
    end
    return nothing
end

function writechains(path::AbstractString, proteinchains::AbstractVector{<:ProteinChain}, format)
    proteinstruc = ProteinStructure(basename(path), proteinchains)
    writechains(path, proteinstruc, format)
end

writechains(path::AbstractString, proteinchain::ProteinChain, format) = writechains(path, [proteinchain], format)

writechains(path::AbstractString, arg) = writechains(path, arg, get_format(path))
