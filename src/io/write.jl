readchains(path::AbstractString, format::Type{<:ProteinFileFormat}) = ProteinStructure(read(path, format))
readchains(path::AbstractString) = readchains(path, get_format(path))

readpdb(path::AbstractString) = readchains(path, PDBFormat)
readcif(path::AbstractString) = readchains(path, MMCIFFormat)

function writepdb(path::AbstractString, chains::AbstractVector{ProteinChain{T}}) where T
    atom_records = BioStructures.AtomRecord[]
    index = 0
    residue_index_global = 0
    for chain in chains
        numbering = hasproperty(chain, :numbering) ? chain.numbering : collect(1:countresidues(chain))
        for residue_index in 1:countresidues(chain)
            resname = get(threeletter_aa_names, chain.sequence[residue_index], "XXX") # threletter_aa_names in residue.jl
            residue_index_global += 1
            residue_atoms = [
                [Atom(atom_name, atom_symbol, coords) for (atom_name, atom_symbol, coords) in zip(BACKBONE_ATOM_NAMES, BACKBONE_ATOM_SYMBOLS, eachcol(view(chain.backbone, :, :, residue_index)))];
                chain.atoms[residue_index]
            ]
            for atom in residue_atoms
                index += 1
                push!(atom_records, BioStructures.AtomRecord(false, index, decode_atom_name(atom.atom_name), ' ', resname, chain.id,
                    numbering[residue_index], ' ', [atom.x, atom.y, atom.z], 1.0, 0.0, atomic_number_to_element_symbol(atom.atomic_number), ""))
            end
        end
    end
    pdblines = BioStructures.pdbline.(atom_records)
    open(path, "w") do io
        for line in pdblines
            println(io, line)
        end
    end
    return nothing
end

writepdb(path::AbstractString, chain::ProteinChain) = writepdb(path, [chain])

writechains(path::AbstractString, chains::AbstractVector{<:ProteinChain}, ::Type{PDBFormat}) = writepdb(path, chains)

function writechains(path::AbstractString, chains::AbstractVector{<:ProteinChain}, format::Type{<:ProteinFileFormat})
    struc = mktempdir() do temp_dir
        temp_path = joinpath(temp_dir, "temp.pdb")
        writechains(temp_path, chains, PDBFormat)
        read(temp_path, PDBFormat) # loads BioStructures.MolecularStructure
    end
    write_function = format == PDBFormat ? BioStructures.writepdb : BioStructures.writemmcif
    write_function(path, struc)
end

writechains(path::AbstractString, chains::AbstractVector{<:ProteinChain}) = writechains(path, chains, get_format(path))

writechains(path, chain::ProteinChain, args...) = writechains(path, [chain], args...)
