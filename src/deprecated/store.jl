import HDF5

struct Lazy{T}
    group::HDF5.Group
end

Base.propertynames(lazy::Lazy, private::Bool=false) = ((private ? (:group,) : ())..., Symbol.(keys(lazy.group))...)

Base.getproperty(lazy::Lazy, property::Symbol) =
    property == :group ? getfield(lazy, property) : readproperty(lazy, Val(property))

Base.setproperty!(lazy::Lazy, property::Symbol, value) =
    property == :group ? setfield!(lazy, property, value) : writeproperty(lazy, Val(property), value)

Base.delete!(lazy::Lazy, property::Symbol) =
    property == :group ? error("group is immutable") : deleteproperty(lazy, Val(property))

include("utils.jl")
include("io.jl")

mutable struct ProteinStructureStore <: AbstractDict{String,ProteinStructure}
    filename::String
    file::HDF5.File
    keys::Set{String}
    mode::String
end

function ProteinStructureStore(filename::AbstractString, mode::AbstractString="cw")
    file = HDF5.h5open(filename, mode)
    store = ProteinStructureStore(filename, file, Set(keys(file)), mode)
    finalizer(close, store)
    return store
end

Base.close(store::ProteinStructureStore) = close(store.file)

Base.keys(store::ProteinStructureStore) = store.keys
Base.length(store::ProteinStructureStore) = length(store.keys)
Base.getindex(store::ProteinStructureStore, key::AbstractString) = read(store.file[key], ProteinStructure)
Base.get(store::ProteinStructureStore, key, default) = key in keys(store) ? store[key] : default

function Base.delete!(store::ProteinStructureStore, key::AbstractString)
    store.mode == "r" && error("$ProteinStructureStore object is read-only.")
    if haskey(store, key)
        HDF5.delete_object(store.file[key])
        delete!(store.keys, key)
    end
    return store
end

function Base.setindex!(store::ProteinStructureStore, value::ProteinStructure, key::AbstractString)
    store.mode == "r" && error("$ProteinStructureStore object is read-only.")
    haskey(store, key) && delete!(store, key)
    group = HDF5.create_group(store.file, key)
    write(group, value)
    push!(store.keys, key)
    return store
end

function Base.iterate(store::ProteinStructureStore, state...)
    iterate_step((key,state)::Tuple{String,Int}) = (key => store[key], state)
    iterate_step(::Nothing) = nothing
    return iterate_step(iterate(keys(store), state...))
end

Base.show(io::IO, store::ProteinStructureStore) = print(io, "$ProteinStructureStore(\"$(store.filename)\", $(store.mode))")
Base.show(io::IO, ::MIME"text/plain", store::ProteinStructureStore) = print(io, summary(store))

Base.open(::Type{ProteinStructureStore}, filename::AbstractString, args...) = ProteinStructureStore(filename, args...)

ProteinStructureStore(f::Function, args...) = open(f, ProteinStructureStore, args...)

function serialize(filename::AbstractString, structures::AbstractVector{<:ProteinStructure})
    ProteinStructureStore(filename, "cw") do store
        for structure in structures
            store[structure.name] = structure
        end
    end
    return filename
end

deserialize(filename::AbstractString) = ProteinStructureStore(collect ∘ values, filename, "r")

function Base.view(store::ProteinStructureStore, name::AbstractString)
    T = eval(Symbol(read(HDF5.attributes(store.file[name])["T"])))
    return Lazy{ProteinStructure{T}}(store.file[name])
end

Base.length(structure::Lazy{<:ProteinStructure}) = read(HDF5.attributes(structure.group)["n_chains"])

Base.getindex(structure::Lazy{ProteinStructure{T}}, index::Integer) where T =
    Lazy{ProteinChain{T}}(structure.group["chains"][string(index)])

Base.length(chain::Lazy{<:ProteinChain}) = read(HDF5.attributes(chain.group)["n_residues"])

# v0.5 compat
function Base.getproperty(lazy::Lazy{T}, property::Symbol) where T<:Union{ProteinStructure,ProteinChain}
    property == :group && return getfield(lazy, :group)
    if haskey(lazy.group, "properties")
        if property in fieldnames(T)
            return invoke(getproperty, Tuple{Lazy,Symbol}, lazy, property)
        else
            return getproperty(readproperty(lazy, Val(:properties)), property)
        end
    else
        return invoke(getproperty, Tuple{Lazy,Symbol}, lazy, property)
    end
end
