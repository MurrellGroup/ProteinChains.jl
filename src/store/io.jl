writeproperty(group::HDF5.Group, name::Symbol, value) = writeproperty(group, Val(name), value)
readproperty(group::HDF5.Group, name::Symbol) = readproperty(group, Val(name))

writeproperty(group::HDF5.Group, ::Val{name}, value) where name = write(group, string(name), value)
readproperty(group::HDF5.Group, ::Val{name}) where name = read(group, string(name))

function writeproperty(group::HDF5.Group, ::Val{:atoms}, atoms::Vector{Vector{Atom{T}}}) where T
    atoms_group = HDF5.create_group(group, "atoms")
    write(atoms_group, "atom_chunk_sizes", map(UInt8∘length, atoms))
    write(atoms_group, "atoms_flattened", reduce(vcat, atoms; init=Atom{T}[]))
    return atoms_group
end

function readproperty(group::HDF5.Group, ::Val{:atoms})
    T = eval(Symbol(read(HDF5.attributes(group)["T"])))
    atoms_group = group["atoms"]
    atom_chunk_sizes = read(atoms_group["atom_chunk_sizes"])
    atoms_flattened = Atom{T}[Atom{T}(atom.name, atom.number, atom.x, atom.y, atom.z) for atom in read(atoms_group["atoms_flattened"])]
    return [atoms_flattened[i-k+1:i] for (i, k) in zip(Iterators.accumulate(+, atom_chunk_sizes; init=0), atom_chunk_sizes)]
end

writeproperty(group::HDF5.Group, ::Val{:numbering}, numbering::Vector{Int32}) = write(group, "numbering", numbers_to_ranges(numbering))
readproperty(group::HDF5.Group, ::Val{:numbering}) = mapreduce(range -> range.start:range.stop, vcat, read(group["numbering"]); init=Int32[])

function writeproperty(group::HDF5.Group, ::Val{:properties}, properties::NamedTuple)
    properties_group = HDF5.create_group(group, "properties")
    subgroups = Dict{Symbol,HDF5.Group}()
    for (name, property) in pairs(properties)
        constructor = typeof(property).name.name
        subgroup = get!(subgroups, constructor) do
            HDF5.create_group(properties_group, string(constructor))
        end
        writeproperty(subgroup, name, unpack(property))
    end
    return properties_group
end

function readproperty(group::HDF5.Group, ::Val{:properties})
    properties_group = group["properties"]
    constructors = Dict{String,Type}(constructor => eval(Symbol(constructor)) for constructor in keys(properties_group))
    return merge(
        (NamedTuple((Symbol(key) => constructors[constructor](readproperty(subgroup, Symbol(key))) for key in keys(subgroup)))
        for (constructor, subgroup) in pairs(properties_group))...
    )
end

function Base.write(group::HDF5.Group, chain::ProteinChain{T}) where T
    HDF5.attributes(group)["T"] = string(T)
    for fieldname in fieldnames(ProteinChain)
        writeproperty(group, fieldname, getproperty(chain, fieldname))
    end
    return group
end

function Base.read(group::HDF5.Group, ::Type{ProteinChain})
    return splat(ProteinChain)(readproperty(group, fieldname) for fieldname in fieldnames(ProteinChain))
end

function Base.write(group::HDF5.Group, structure::ProteinStructure{T}) where T
    HDF5.attributes(group)["T"] = string(T)
    HDF5.attributes(group)["n_residues"] = sum(length, structure)
    HDF5.attributes(group)["n_chains"] = length(structure)

    group["name"] = structure.name
    group["atoms"] = structure.atoms

    chains_group = HDF5.create_group(group, "chains")

    for (i, chain) in enumerate(structure.chains)
        chain_group = HDF5.create_group(chains_group, string(i))
        write(chain_group, chain)
    end

    return group
end

function Base.read(group::HDF5.Group, ::Type{ProteinStructure})
    T = eval(Symbol(read(HDF5.attributes(group)["T"])))
    name = read(group["name"])
    atoms = [Atom(atom.name, atom.number, atom.x, atom.y, atom.z) for atom in read(group["atoms"])]

    chains_group = group["chains"]
    chains = ProteinChain{T}[]
    for key in keys(chains_group)
        chain = read(chains_group[key], ProteinChain)
        push!(chains, chain)
    end

    return ProteinStructure(name, atoms, chains)
end
