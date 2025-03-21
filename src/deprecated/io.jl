# v0.5 compat

const _NEW_CONSTRUCTOR_MAP = Dict("StandardProperty" => :identity, "IndexableProperty" => :Indexable)

function readproperty(lazy::Lazy{T}, ::Val{:properties}) where T
    properties_group = lazy.group["properties"]
    constructors = Dict(constructor => get(_NEW_CONSTRUCTOR_MAP, constructor, Symbol(constructor)) for constructor in keys(properties_group))
    isempty(constructors) && return (;)
    return merge(
        (NamedTuple((Symbol(key) => eval(constructors[constructor])(readproperty(Lazy{T}(subgroup), Val(Symbol(key)))) for key in keys(subgroup)))
        for (constructor, subgroup) in pairs(properties_group))...
    )
end

# generic property read/write

function writeproperty(lazy::Lazy, ::Val{name}, value) where name
    write(lazy.group, string(name), value)
end

function readproperty(lazy::Lazy, ::Val{name}) where name
    str = string(name)
    value = read(lazy.group, str)
    attr = HDF5.attributes(lazy.group)
    f = haskey(attr, str) ? eval(Symbol(read(attr[str]))) : identity
    return f(value)
end

function deleteproperty(lazy::Lazy, ::Val{name}) where name
    if haskey(lazy.group, string(name))
        HDF5.delete_object(lazy.group, string(name))
        dict = HDF5.AttributeDict(lazy.group)
        haskey(dict, string(name)) && delete!(dict, string(name))
    end
end

function writeproperty(lazy::Lazy, ::Val{name}, value::AbstractString) where name
    HDF5.attributes(lazy.group)[string(name)] = "String"
    writeproperty(lazy, Val(name), collect(codeunits(value)))
end

function writeproperty(lazy::Lazy, ::Val{name}, property::AbstractProperty) where name
    HDF5.attributes(lazy.group)[string(name)] = string(typeof(property))
    write(lazy.group, string(name), unwrap(property))
end

import JSON
const _json = JSON.parse ∘ String
function writeproperty(lazy::Lazy, ::Val{name}, x::Union{AbstractDict, AbstractVector{<:Union{AbstractDict, AbstractVector}}}) where name
    HDF5.attributes(lazy.group)[string(name)] = "_json"
    write(lazy.group, string(name), collect(codeunits(JSON.json(x))))
end

## chain

function writeproperty(lazy::Lazy{ProteinChain{T}}, ::Val{:atoms}, atoms::Vector{Vector{Atom{T}}}) where T
    atoms_group = HDF5.create_group(lazy.group, "atoms")
    write(atoms_group, "atom_chunk_sizes", map(UInt8∘length, atoms))
    write(atoms_group, "atoms_flattened", reduce(vcat, atoms; init=Atom{T}[]))
end

new_name(name::UInt32) = convert(ProteinChains.AtomName, String(reinterpret(UInt8, [name])))

function readproperty(lazy::Lazy{ProteinChain{T}}, ::Val{:atoms}) where T
    atoms_group = lazy.group["atoms"]
    atom_chunk_sizes = read(atoms_group["atom_chunk_sizes"])
    atoms_flattened = map(read(atoms_group["atoms_flattened"])) do atom
        Atom(new_name(atom.name), atom.number, atom.x, atom.y, atom.z)
    end
    return [atoms_flattened[i-k+1:i] for (i, k) in zip(Iterators.accumulate(+, atom_chunk_sizes; init=0), atom_chunk_sizes)]
end

writeproperty(lazy::Lazy{ProteinChain{T}}, ::Val{:ins_codes}, ins_codes::String) where T =
    write(lazy.group, "ins_codes", encode_ins_codes(ins_codes))
readproperty(lazy::Lazy{ProteinChain{T}}, ::Val{:ins_codes}) where T =
    decode_ins_codes(read(lazy.group["ins_codes"]), read(HDF5.attributes(lazy.group)["n_residues"]))

writeproperty(lazy::Lazy{ProteinChain{T}}, ::Val{:numbering}, numbering::Vector{<:Integer}) where T =
    write(lazy.group, "numbering", numbers_to_ranges(convert(Vector{Int32}, numbering)))
readproperty(lazy::Lazy{ProteinChain{T}}, ::Val{:numbering}) where T =
    mapreduce(range -> range.start:range.stop, vcat, read(lazy.group["numbering"]); init=Int[])

function Base.write(group::HDF5.Group, chain::ProteinChain{T}) where T
    HDF5.attributes(group)["n_residues"] = length(chain)
    for name in propertynames(chain)
        writeproperty(Lazy{ProteinChain{T}}(group), Val(name), getproperty(chain, name))
    end
    return group
end

function Base.read(group::HDF5.Group, ::Type{ProteinChain{T}}) where T
    fields = fieldnames(ProteinChain)
    properties = Symbol.(keys(group))
    args = Iterators.map(Iterators.filter(in(properties), fields)) do name
        readproperty(Lazy{ProteinChain{T}}(group), Val(name))
    end
    kws = if haskey(group, "properties")
        # v0.5 compat
        merge(readproperty(Lazy{ProteinChain{T}}(group), Val(:properties)), (; ins_codes=Indexable(collect(readproperty(Lazy{ProteinChain{T}}(group), Val(:ins_codes))))))
    else
        Iterators.map(setdiff(properties, fields)) do name
            name => readproperty(Lazy{ProteinChain{T}}(group), Val(name))
        end
    end
    return ProteinChain(args...; kws...)
end

## structure

function readproperty(lazy::Lazy{ProteinStructure{T}}, ::Val{:atoms}) where T
    atoms = read(lazy.group["atoms"])
    return collect(reinterpret(Atom{T}, atoms))
end

function writeproperty(lazy::Lazy{ProteinStructure{T}}, ::Val{:chains}, chains::Vector{<:ProteinChain}) where T
    chains_group = HDF5.create_group(lazy.group, "chains")
    for (i, chain) in enumerate(chains)
        chain_group = HDF5.create_group(chains_group, string(i))
        write(chain_group, chain)
    end
end

function readproperty(lazy::Lazy{ProteinStructure{T}}, ::Val{:chains}) where T
    chains_group = lazy.group["chains"]
    chain_numbers = sort(collect(keys(chains_group)), by=i->parse(Int, i))
    return [read(chains_group[key], ProteinChain{T}) for key in chain_numbers]
end

function Base.write(group::HDF5.Group, structure::ProteinStructure{T}) where T
    HDF5.attributes(group)["T"] = string(T)
    HDF5.attributes(group)["n_chains"] = length(structure)
    for name in propertynames(structure)
        writeproperty(Lazy{ProteinStructure{T}}(group), Val(name), getproperty(structure, name))
    end
    return group
end

function Base.read(group::HDF5.Group, ::Type{ProteinStructure{T}}) where T
    fields = fieldnames(ProteinStructure)
    properties = Symbol.(keys(group))
    args = Iterators.map(Iterators.filter(in(properties), fields)) do name
        readproperty(Lazy{ProteinStructure{T}}(group), Val(name))
    end
    kws = if haskey(group, "properties")
        # v0.5 compat
        #@warn "Reading deprecated ProteinStructureStore file format version. Please reserialize the data to a new file."
        readproperty(Lazy{ProteinStructure{T}}(group), Val(:properties))
    else
        Iterators.map(setdiff(properties, fields)) do name
            name => readproperty(Lazy{ProteinStructure{T}}(group), Val(name))
        end
    end
    return ProteinStructure(args...; kws...)
end

Base.read(group::HDF5.Group, ::Type{ProteinStructure}) = read(group, ProteinStructure{eval(Symbol(read(HDF5.attributes(group)["T"])))})

Base.read(lazy::Lazy{T}) where T = read(lazy.group, T)
