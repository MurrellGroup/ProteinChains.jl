### FOR BACKWARD READING COMPATIBILITY

const _NEW_CONSTRUCTOR_MAP = Dict("StandardProperty" => :identity, "IndexableProperty" => :Indexable)

function readproperty(group::HDF5.Group, T::Type, ::Val{:properties})
    properties_group = group["properties"]
    constructors = Dict(constructor => get(_NEW_CONSTRUCTOR_MAP, constructor, Symbol(constructor)) for constructor in keys(properties_group))
    isempty(constructors) && return (;)
    return merge(
        (NamedTuple((Symbol(key) => eval(constructors[constructor])(readproperty(subgroup, T, Val(Symbol(key)))) for key in keys(subgroup)))
        for (constructor, subgroup) in pairs(properties_group))...
    )
end


function writeproperty(group::HDF5.Group, ::Type, ::Val{name}, value) where name
    haskey(group, string(name)) && HDF5.delete_object(group[string(name)])
    return write(group, string(name), value)
end

function writeproperty(group::HDF5.Group, T::Type, ::Val{name}, property::AbstractProperty) where name
    HDF5.attributes(group)[string(name)] = string(typeof(property))
    return writeproperty(group, T, Val(name), unwrap(property))
end

function readproperty(group::HDF5.Group, ::Type, ::Val{name}) where name
    str = string(name)
    value = read(group, str)
    attr = HDF5.attributes(group)
    f = haskey(attr, str) ? eval(Symbol(read(attr[str]))) : identity
    return f(value)
end

import JSON
const _json = JSON.parse
function writeproperty(group::HDF5.Group, ::Type, ::Val{name}, dict::AbstractDict) where name
    HDF5.attributes(group)[string(name)] = "_json"
    write(group, string(name), JSON.json(dict))
end

## chain

function writeproperty(group::HDF5.Group, ::Type{ProteinChain{T}}, ::Val{:atoms}, atoms::Vector{Vector{Atom{T}}}) where T
    atoms_group = HDF5.create_group(group, "atoms")
    write(atoms_group, "atom_chunk_sizes", map(UInt8∘length, atoms))
    write(atoms_group, "atoms_flattened", reduce(vcat, atoms; init=Atom{T}[]))
    return atoms_group
end

function readproperty(group::HDF5.Group, ::Type{ProteinChain{T}}, ::Val{:atoms}) where T
    atoms_group = group["atoms"]
    atom_chunk_sizes = read(atoms_group["atom_chunk_sizes"])
    atoms_flattened = collect(reinterpret(Atom{T}, read(atoms_group["atoms_flattened"])))
    return [atoms_flattened[i-k+1:i] for (i, k) in zip(Iterators.accumulate(+, atom_chunk_sizes; init=0), atom_chunk_sizes)]
end

writeproperty(group::HDF5.Group, ::Type{ProteinChain{T}}, ::Val{:ins_codes}, ins_codes::String) where T =
    write(group, "ins_codes", encode_ins_codes(ins_codes))
readproperty(group::HDF5.Group, ::Type{ProteinChain{T}}, ::Val{:ins_codes}) where T =
    decode_ins_codes(read(group["ins_codes"]), read(HDF5.attributes(group)["n_residues"]))

writeproperty(group::HDF5.Group, ::Type{ProteinChain{T}}, ::Val{:numbering}, numbering::Vector{<:Integer}) where T =
    write(group, "numbering", numbers_to_ranges(convert(Vector{Int32}, numbering)))
readproperty(group::HDF5.Group, ::Type{ProteinChain{T}}, ::Val{:numbering}) where T =
    mapreduce(range -> range.start:range.stop, vcat, read(group["numbering"]); init=Int[])

function Base.write(group::HDF5.Group, chain::ProteinChain{T}) where T
    HDF5.attributes(group)["n_residues"] = length(chain)
    for name in propertynames(chain)
        writeproperty(group, ProteinChain{T}, Val(name), getproperty(chain, name))
    end
    return group
end

function Base.read(group::HDF5.Group, ::Type{ProteinChain{T}}) where T
    fields = fieldnames(ProteinChain)
    properties = Symbol.(keys(group))
    args = Iterators.map(Iterators.filter(in(properties), fields)) do name
        readproperty(group, ProteinChain{T}, Val(name))
    end
    kwargs = if haskey(group, "properties")
        # for backward compatibility
        readproperty(group, ProteinChain{T}, Val(:properties))
    else
        Iterators.map(setdiff(properties, fields)) do name
            name => readproperty(group, ProteinChain{T}, Val(name))
        end
    end
    return ProteinChain(args...; kwargs...)
end

## structure

writeproperty(group::HDF5.Group, ::Type{ProteinStructure{T}}, ::Val{:atoms}, atoms::Vector{Atom{T}}) where T =
    write(group, "atoms", atoms)

readproperty(group::HDF5.Group, ::Type{ProteinStructure{T}}, ::Val{:atoms}) where T =
   collect(reinterpret(Atom{T}, read(group["atoms"])))

function writeproperty(group::HDF5.Group, ::Type{<:ProteinStructure}, ::Val{:chains}, chains::Vector{<:ProteinChain})
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
    for name in propertynames(structure)
        writeproperty(group, ProteinStructure{T}, Val(name), getproperty(structure, name))
    end
    return group
end

function Base.read(group::HDF5.Group, ::Type{ProteinStructure{T}}) where T
    fields = fieldnames(ProteinStructure)
    properties = Symbol.(keys(group))
    args = Iterators.map(Iterators.filter(in(properties), fields)) do name
        readproperty(group, ProteinStructure{T}, Val(name))
    end
    kwargs = if haskey(group, "properties")
        # for backward compatibility
        @warn "Reading deprecated ProteinStructureStore file format. Please reserialize the data to a new file."
        readproperty(group, ProteinStructure{T}, Val(:properties))
    else
        Iterators.map(setdiff(properties, fields)) do name
            name => readproperty(group, ProteinStructure{T}, Val(name))
        end
    end
    return ProteinStructure(args...; kwargs...)
end

Base.read(group::HDF5.Group, ::Type{ProteinStructure}) = read(group, ProteinStructure{eval(Symbol(read(HDF5.attributes(group)["T"])))})
