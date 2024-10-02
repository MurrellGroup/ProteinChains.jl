import HDF5

function Base.write(parent::Union{HDF5.File,HDF5.Group}, chain::ProteinChain{T}, path::AbstractString=chain.id) where T
    chain_group = HDF5.create_group(parent, path)
    HDF5.attributes(chain_group)["T"] = string(T)

    chain_group["id"] = chain.id
    chain_group["atom_chunk_sizes"] = map(UInt8∘length, chain.atoms)
    chain_group["atoms_flattened"] = reduce(vcat, chain.atoms)
    chain_group["sequence"] = chain.sequence
    chain_group["numbering"] = map(UInt32, chain.numbering)

    indexable = HDF5.create_group(chain_group, "indexable")
    persistent = HDF5.create_group(chain_group, "persistent")
    for (name, property) in pairs(chain.properties)
        g = property isa ResidueProperty ? indexable : persistent
        g[string(name)] = property[]
    end

    return parent
end

function Base.read(group::Union{Union{HDF5.File,HDF5.Group}}, ::Type{ProteinChain})
    T = eval(Symbol(read(HDF5.attributes(group)["T"])))
    id = read(group["id"])
    atom_chunk_sizes = read(group["atom_chunk_sizes"])
    atoms_flattened = [Atom(nt.name, nt.number, nt.x, nt.y, nt.z) for nt in read(group["atoms_flattened"])]
    sequence = read(group["sequence"])
    numbering = read(group["numbering"])

    atoms = Vector{Atom{T}}[]
    for (i, k) in zip(Iterators.accumulate(+, atom_chunk_sizes), atom_chunk_sizes)
        residue_atoms = atoms_flattened[i-k+1:i]
        push!(atoms, residue_atoms)
    end

    indexable = group["indexable"]
    persistent = group["persistent"]
    properties = merge(
        NamedTuple((Symbol(key) => ResidueProperty(read(indexable[key])) for key in keys(group["indexable"]))),
        NamedTuple((Symbol(key) => ChainProperty(read(persistent[key])) for key in keys(group["persistent"])))
    )

    ProteinChain(id, atoms, sequence, numbering, properties)
end

function Base.write(parent::Union{HDF5.File,HDF5.Group}, structure::ProteinStructure{T}, path::AbstractString=structure.name) where T
    structure_group = HDF5.create_group(parent, path)
    HDF5.attributes(structure_group)["T"] = string(T)

    structure_group["name"] = structure.name
    structure_group["atoms"] = structure.atoms
    structure_group["numbering"] = structure.numbering

    chains_group = HDF5.create_group(structure_group, "chains")

    for (chain, number) in zip(structure.chains, structure.numbering)
        write(chains_group, chain, string(number))
    end

    return parent
end

function Base.read(group::Union{Union{HDF5.File,HDF5.Group}}, ::Type{ProteinStructure})
    T = eval(Symbol(read(HDF5.attributes(group)["T"])))
    name = read(group["name"])
    atoms = [Atom(nt.name, nt.number, nt.x, nt.y, nt.z) for nt in read(group["atoms"])]
    numbering = read(group["numbering"])

    chains_group = group["chains"]
    chains = ProteinChain{T}[]
    for i in numbering
        chain = read(chains_group[string(i)], ProteinChain)
        push!(chains, chain)
    end

    return ProteinStructure(name, atoms, chains, numbering)
end

function serialize(path::AbstractString, structures::AbstractVector{<:ProteinStructure}; mode="w")
    HDF5.h5open(path, mode) do f
        foreach(structure -> write(f, structure), structures)
    end
    return path
end

function deserialize(path::AbstractString)
    structures = ProteinStructure[]
    HDF5.h5open(path, "r") do f
        for key in keys(f)
            push!(structures, read(f[key], ProteinStructure))
        end
    end
    return structures
end
