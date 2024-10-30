function writeproperty(group::HDF5.Group, ::Type, ::Val{name}, value) where name
    haskey(group, string(name)) && HDF5.delete_object(group[string(name)])
    return write(group, string(name), value)
end

readproperty(group::HDF5.Group, ::Type, ::Val{name}) where name = read(group, string(name))

function writeproperty(group::HDF5.Group, T::Type, ::Val{:properties}, properties::NamedTuple)
    properties = namedproperties(properties)
    properties_group = haskey(group, "properties") ? group["properties"] : HDF5.create_group(group, "properties")
    subgroups = Dict{Symbol,HDF5.Group}(Symbol(constructor) => subgroup for (constructor, subgroup) in pairs(properties_group))
    for (name, property) in pairs(properties)
        constructor = typeof(property).name.name
        subgroup = get!(subgroups, constructor) do
            HDF5.create_group(properties_group, string(constructor))
        end
        writeproperty(subgroup, T, Val(name), unpack(property))
    end
    return properties_group
end

function readproperty(group::HDF5.Group, T::Type, ::Val{:properties})
    properties_group = group["properties"]
    constructors = Dict{String,Type}(constructor => eval(Symbol(constructor)) for constructor in keys(properties_group))
    isempty(constructors) && return (;)
    return merge(
        (NamedTuple((Symbol(key) => constructors[constructor](readproperty(subgroup, T, Val(Symbol(key)))) for key in keys(subgroup)))
        for (constructor, subgroup) in pairs(properties_group))...
    )
end

function deleteproperty(group::HDF5.Group, ::Type, ::Val{:properties}, names::Symbol...)
    properties_group = group["properties"]
    for name in names
        for subgroup in values(properties_group)
            haskey(subgroup, string(name)) && HDF5.delete_object(subgroup[string(name)])
        end
    end
    return group
end

## chain properties

function writeproperty(group::HDF5.Group, ::Type{ProteinChain{T}}, ::Val{:atoms}, atoms::Vector{Vector{Atom{T}}}) where T
    atoms_group = HDF5.create_group(group, "atoms")
    write(atoms_group, "atom_chunk_sizes", map(UInt8∘length, atoms))
    write(atoms_group, "atoms_flattened", reduce(vcat, atoms; init=Atom{T}[]))
    return atoms_group
end

function readproperty(group::HDF5.Group, ::Type{ProteinChain{T}}, ::Val{:atoms}) where T
    atoms_group = group["atoms"]
    atom_chunk_sizes = read(atoms_group["atom_chunk_sizes"])
    atoms_flattened = [Atom{T}(atom.name, atom.number, atom.x, atom.y, atom.z) for atom in read(atoms_group["atoms_flattened"])]
    return [atoms_flattened[i-k+1:i] for (i, k) in zip(Iterators.accumulate(+, atom_chunk_sizes; init=0), atom_chunk_sizes)]
end

writeproperty(group::HDF5.Group, ::Type{ProteinChain{T}}, ::Val{:ins_codes}, ins_codes::String) where T =
    write(group, "ins_codes", encode_ins_codes(ins_codes))
readproperty(group::HDF5.Group, ::Type{ProteinChain{T}}, ::Val{:ins_codes}) where T =
    decode_ins_codes(read(group["ins_codes"]), read(HDF5.attributes(group)["n_residues"]))

writeproperty(group::HDF5.Group, ::Type{ProteinChain{T}}, ::Val{:numbering}, numbering::Vector{Int32}) where T =
    write(group, "numbering", numbers_to_ranges(numbering))
readproperty(group::HDF5.Group, ::Type{ProteinChain{T}}, ::Val{:numbering}) where T =
    mapreduce(range -> range.start:range.stop, vcat, read(group["numbering"]); init=Int32[])

function Base.write(group::HDF5.Group, chain::ProteinChain{T}) where T
    HDF5.attributes(group)["n_residues"] = length(chain)
    for fieldname in fieldnames(ProteinChain)
        writeproperty(group, ProteinChain{T}, Val(fieldname), getproperty(chain, fieldname))
    end
    return group
end

Base.read(group::HDF5.Group, ::Type{ProteinChain{T}}) where T =
    splat(ProteinChain{T})(readproperty(group, ProteinChain{T}, Val(fieldname)) for fieldname in fieldnames(ProteinChain))

## structure properties

writeproperty(group::HDF5.Group, ::Type{ProteinStructure{T}}, ::Val{:atoms}, atoms::Vector{Atom{T}}) where T =
    write(group, "atoms", atoms)

readproperty(group::HDF5.Group, ::Type{ProteinStructure{T}}, ::Val{:atoms}) where T =
    [Atom{T}(atom.name, atom.number, atom.x, atom.y, atom.z) for atom in read(group["atoms"])]

function writeproperty(group::HDF5.Group, ::Type{ProteinStructure{T}}, ::Val{:chains}, chains::Vector{ProteinChain{T}}) where T
    chains_group = HDF5.create_group(group, "chains")
    for (i, chain) in enumerate(chains)
        chain_group = HDF5.create_group(chains_group, string(i))
        write(chain_group, chain)
    end
    return group
end

function readproperty(group::HDF5.Group, ::Type{ProteinStructure{T}}, ::Val{:chains}) where T
    chains_group = group["chains"]
    chain_numbers = sort(collect(keys(chains_group)), by=i->parse(Int, i))
    return [read(chains_group[key], ProteinChain{T}) for key in chain_numbers]
end

function Base.write(group::HDF5.Group, structure::ProteinStructure{T}) where T
    HDF5.attributes(group)["T"] = string(T)
    HDF5.attributes(group)["n_residues"] = sum(length, structure)
    HDF5.attributes(group)["n_chains"] = length(structure)
    for fieldname in fieldnames(ProteinStructure)
        writeproperty(group, ProteinStructure{T}, Val(fieldname), getproperty(structure, fieldname))
    end
    return group
end

Base.read(group::HDF5.Group, ::Type{ProteinStructure{T}}) where T =
    splat(ProteinStructure{T})(readproperty(group, ProteinStructure{T}, Val(fieldname)) for fieldname in fieldnames(ProteinStructure))

Base.read(group::HDF5.Group, ::Type{ProteinStructure}) = read(group, ProteinStructure{eval(Symbol(read(HDF5.attributes(group)["T"])))})
