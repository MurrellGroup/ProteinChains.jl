using LinearAlgebra: cross, normalize

_normalize(A::AbstractArray; dims=1) = mapslices(normalize, A; dims)

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

function get_oxygen_positions(backbone::Array{T,3}) where T<:Real
    Ca_pos, C_pos, N_pos = backbone[:, 2, 1:end-1], backbone[:, 3, 1:end-1], backbone[:, 1, 2:end]
    CaC_vec = C_pos - Ca_pos
    NC_vec = C_pos - N_pos
    CO_vec = T(1.23) * _normalize(_normalize(CaC_vec) + _normalize(NC_vec))
    O_pos = C_pos + CO_vec
    return [O_pos estimate_last_oxygen_position(backbone)]
end

function BioStructures.Chain(proteinchain::ProteinChain, model::BioStructures.Model)
    numbering = proteinchain.numbering
    residue_list = Vector{String}()
    residues = Dict{String, BioStructures.AbstractResidue}()
    chain = BioStructures.Chain(proteinchain.id, residue_list, residues, model)

    # some secondary structure visualization tools require oxygen atoms to be present,
    # so these get used if a residue is missing the oxygen atom
    oxygen_name = encode_atom_name("O", "O")
    backbone_coords = get_backbone(proteinchain)
    oxygen_coords = get_oxygen_positions(backbone_coords)

    atom_serial = 0
    for residue_index in 1:length(proteinchain)
        atom_list = Vector{String}()
        atoms = Dict{String, BioStructures.AbstractAtom}()
        resname = threeletter_resname(proteinchain.sequence[residue_index])
        number = numbering[residue_index]
        residue = BioStructures.Residue(resname, number, ' ', false, atom_list, atoms, chain, '-') # TODO: secondary structure

        # each residue has at least N, CA, C atoms
        for atom in proteinchain.atoms[residue_index]
            atom_serial += 1
            name = decode_atom_name(atom.name)
            coords = [atom.x, atom.y, atom.z]
            element = number_to_element_symbol(atom.number)
            atom = BioStructures.Atom(atom_serial, name, ' ', coords, 1.0, 0.0, element, "  ", residue)
            push!(atom_list, atom.name)
            atoms[atom.name] = atom
        end

        if !any(atom -> atom.name == oxygen_name, proteinchain.atoms[residue_index])
            atom_serial += 1
            coords = oxygen_coords[:, residue_index]
            atom = BioStructures.Atom(atom_serial, "O", ' ', coords, 1.0, 0.0, "O", "  ", residue)
            insert!(atom_list, 4, atom.name)
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
        modelnum = 1 # only supporting one model for now
        model = get!(models, modelnum) do
            BioStructures.Model(modelnum, Dict{String, BioStructures.Chain}(), struc)
        end
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
end

function Base.write(path::AbstractString, proteinchains::AbstractVector{<:ProteinChain{T}}, format=get_format(path)) where T
    proteinstruc = ProteinStructure(basename(path), Atom{T}[], proteinchains)
    write(path, proteinstruc, format)
end

Base.write(path::AbstractString, proteinchain::ProteinChain, format=get_format(path)) = write(path, [proteinchain], format)

writecif(path::AbstractString, x) = write(path, x, MMCIFFormat)
writepdb(path::AbstractString, x) = write(path, x, PDBFormat)
