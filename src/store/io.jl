include("utils.jl")

function writeh5(group::HDF5.Group, chain::ProteinChain{T}) where T
    HDF5.attributes(group)["T"] = string(T)

    group["id"] = chain.id
    group["atom_chunk_sizes"] = map(UInt8∘length, chain.atoms)
    group["atoms_flattened"] = reduce(vcat, chain.atoms; init=Atom{T}[])
    group["sequence"] = chain.sequence
    group["numbering"] = numbers_to_ranges(map(Int32, chain.numbering))

    indexable = HDF5.create_group(group, "indexable")
    persistent = HDF5.create_group(group, "persistent")
    for (name, property) in pairs(chain.properties)
        g = property isa IndexableProperty ? indexable : persistent
        g[string(name)] = property[]
    end

    return group
end

function readh5(group::HDF5.Group, ::Type{ProteinChain})
    T = eval(Symbol(read(HDF5.attributes(group)["T"])))
    id = read(group["id"])
    atom_chunk_sizes = read(group["atom_chunk_sizes"])
    atoms_flattened = [Atom(atom.name, atom.number, atom.x, atom.y, atom.z) for atom in read(group["atoms_flattened"])]
    sequence = read(group["sequence"])
    numbering = mapreduce(range -> range.start:range.stop, vcat, read(group["numbering"]); init=Int32[])

    atoms = Vector{Atom{T}}[]
    for (i, k) in zip(Iterators.accumulate(+, atom_chunk_sizes; init=0), atom_chunk_sizes)
        push!(atoms, atoms_flattened[i-k+1:i])
    end

    indexable = group["indexable"]
    persistent = group["persistent"]
    properties = merge(
        NamedTuple((Symbol(key) => IndexableProperty(read(indexable[key])) for key in keys(group["indexable"]))),
        NamedTuple((Symbol(key) => PersistentProperty(read(persistent[key])) for key in keys(group["persistent"])))
    )

    return ProteinChain(id, atoms, sequence, numbering, properties)
end

function writeh5(group::HDF5.Group, structure::ProteinStructure{T}) where T
    HDF5.attributes(group)["T"] = string(T)

    group["name"] = structure.name
    group["atoms"] = structure.atoms

    chains_group = HDF5.create_group(group, "chains")

    for (i, chain) in enumerate(structure.chains)
        chain_group = HDF5.create_group(chains_group, string(i))
        writeh5(chain_group, chain)
    end

    return group
end

function readh5(group::HDF5.Group, ::Type{ProteinStructure})
    T = eval(Symbol(read(HDF5.attributes(group)["T"])))
    name = read(group["name"])
    atoms = [Atom(atom.name, atom.number, atom.x, atom.y, atom.z) for atom in read(group["atoms"])]

    chains_group = group["chains"]
    chains = ProteinChain{T}[]
    for key in keys(chains_group)
        chain = readh5(chains_group[key], ProteinChain)
        push!(chains, chain)
    end

    return ProteinStructure(name, atoms, chains)
end
