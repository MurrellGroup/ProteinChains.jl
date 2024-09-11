using LinearAlgebra: cross, normalize

using AssigningSecondaryStructure: get_oxygen_positions

function estimate_last_oxygen_position(backbone::Array{T,3}) where T<:Real
    N_pos, Ca_pos, C_pos = eachcol(backbone[:,:,end])
    angle = 2π/3
    bond_length = 1.2
    v = Ca_pos - C_pos
    w = cross(v, Ca_pos - N_pos)
    u = cross(w, v)
    O_pos = C_pos + normalize(u)*bond_length*cos(angle - 0.5π) - normalize(v)*bond_length*sin(angle - 0.5π)
    return O_pos
end

function BioStructures.Chain(proteinchain::AbstractProteinChain, model::BioStructures.Model)
    numbering = proteinchain.numbering
    residue_list = Vector{String}()
    residues = Dict{String, BioStructures.AbstractResidue}()
    chain = BioStructures.Chain(proteinchain.id, residue_list, residues, model)

    # some visualization tools require oxygen atoms to be present, so these get used if a residue is missing the oxygen atom
    oxygen_name = encode_atom_name("O", "O")
    oxygen_coords = get_oxygen_positions(proteinchain.backbone)
    oxygen_coords[:,end] = estimate_last_oxygen_position(proteinchain.backbone)

    atom_serial = 0
    for residue_index in 1:length(proteinchain)
        atom_list = Vector{String}()
        atoms = Dict{String, BioStructures.AbstractAtom}()
        resname = threeletter_resname(proteinchain.sequence[residue_index])
        number = numbering[residue_index]
        residue = BioStructures.Residue(resname, number, ' ', false, atom_list, atoms, chain, '-') # TODO: secondary structure

        for (i, (name, element)) in enumerate(zip(BACKBONE_ATOM_NAMES, BACKBONE_ATOM_SYMBOLS))
            atom_serial += 1
            coords = proteinchain.backbone[:, i, residue_index]
            atom = BioStructures.Atom(atom_serial, name, ' ', coords, 1.0, 0.0, element, "  ", residue)
            push!(atom_list, atom.name)
            atoms[atom.name] = atom
        end

        if !any(atom -> atom.name == oxygen_name, proteinchain.atoms[residue_index])
            atom_serial += 1
            coords = oxygen_coords[:, residue_index]
            atom = BioStructures.Atom(atom_serial, "O", ' ', coords, 1.0, 0.0, "O", "  ", residue)
            push!(atom_list, atom.name)
            atoms[atom.name] = atom
        end

        for atom in proteinchain.atoms[residue_index]
            atom_serial += 1
            name = decode_atom_name(atom.name)
            coords = [atom.x, atom.y, atom.z]
            element = atomic_number_to_element_symbol(atom.atomic_number)
            atom = BioStructures.Atom(atom_serial, name, ' ', coords, 1.0, 0.0, element, "  ", residue)
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

function Base.write(path::AbstractString, proteinstruc::ProteinStructure, format::Type{<:ProteinFileFormat})
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

function Base.write(path::AbstractString, proteinchains::AbstractVector{<:AbstractProteinChain}, format=get_format(path))
    proteinstruc = ProteinStructure(basename(path), proteinchains)
    write(path, proteinstruc, format)
end

Base.write(path::AbstractString, proteinchain::AbstractProteinChain, format=get_format(path)) = write(path, [proteinchain], format)

writecif(path::AbstractString, x) = write(path, x, MMCIFFormat)
writepdb(path::AbstractString, x) = write(path, x, PDBFormat)
