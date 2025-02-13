import HDF5

include("utils.jl")
include("io.jl")

"""
    ProteinStructureStore <: AbstractDict{String,ProteinStructure}

A mutable struct representing a store for protein structures using HDF5 file format.

The struct implements the AbstractDict interface, allowing for dictionary-like operations.

A `ProteinStructureStore` gets closed when there no longer exists a program-accessible reference to it.

## Examples

```jldoctest
julia> store = ProteinStructureStore("store.h5")
ProteinStructureStore with 0 entries

julia> store["3HFM"] = pdb"3HFM"
[ Info: Downloading file from PDB: 3HFM
3-element ProteinStructure "3HFM.cif"
 215-residue ProteinChain{Float64} (H)
 214-residue ProteinChain{Float64} (L)
 129-residue ProteinChain{Float64} (Y)

julia> store
ProteinStructureStore with 1 entry

julia> keys(store)
Set{String} with 1 element:
  "3HFM"

julia> delete!(store, "3HFM")
ProteinStructureStore with 0 entries
```
"""
mutable struct ProteinStructureStore <: AbstractDict{String,ProteinStructure}
    filename::String
    file::HDF5.File
    keys::Set{String}
    mode::String
end

"""
    ProteinStructureStore(filename, mode="cw")

Open or create an HDF5 file as a `ProteinStructureStore` where `mode` is one of:
- "r" read only
- "r+" read and write
- "cw" read and write, create file if not existing, do not truncate
- "w" read and write, create a new file (destroys any existing contents)
"""
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

"""
    ProteinStructureStore(f::Function, filename, mode="cw")
"""
ProteinStructureStore(f::Function, args...) = open(f, ProteinStructureStore, args...)

"""
    serialize(filename::AbstractString, structures::AbstractVector{<:ProteinStructure})

Serialize a vector of `ProteinStructure` objects to an HDF5 file.
This function creates a new `ProteinStructureStore` and writes each structure in the input vector to it.
Each structure is stored using its name as the key.
"""
function serialize(filename::AbstractString, structures::AbstractVector{<:ProteinStructure})
    ProteinStructureStore(filename, "cw") do store
        for structure in structures
            store[structure.name] = structure
        end
    end
    return filename
end

"""
    deserialize(filename::AbstractString)

Deserialize `ProteinStructure` objects from an HDF5 file.
Returns a `Vector{ProteinStructure}` of all structures stored in the file.
"""
deserialize(filename::AbstractString) = ProteinStructureStore(collect ∘ values, filename, "r")

struct Lazy{T}
    group::HDF5.Group
end

Base.length(chain::Lazy{<:ProteinChain}) = read(HDF5.attributes(chain.group)["n_residues"])
Base.length(structure::Lazy{<:ProteinStructure}) = read(HDF5.attributes(structure.group)["n_chains"])

function Base.getproperty(lazy::Lazy{T}, property::Symbol) where T
    property == :group && return getfield(lazy, property)
    haskey(lazy.group, string(property)) && return readproperty(lazy.group, T, Val(property))
    # v0.5 compat
    haskey(lazy.group, "properties") && return getproperty(readproperty(lazy.group, T, Val(:properties)), property)
end

Base.getindex(structure::Lazy{ProteinStructure{T}}, index::Integer) where T = Lazy{ProteinChain{T}}(structure.group["chains"][string(index)])

"""
    Base.view(store::ProteinStructureStore, name::AbstractString)

Return a lazy view of a structure in the store, allowing for partial loading of structures.

```jldoctest
julia> store = ProteinStructureStore("store.pss");

julia> store["3NIR"] = pdb"3NIR"
[ Info: Downloading file from PDB: 3NIR
1-element ProteinStructure "3NIR.cif":
 46-residue ProteinChain{Float64} (A)

julia> view(store, "3NIR").name
"3NIR.cif"

julia> view(store, "3NIR")[1].sequence
"TTCCPSIVARSNFNVCRLPGTPEALCATYTGCIIIPGATCPGDYAN"
```
"""
function Base.view(store::ProteinStructureStore, name::AbstractString)
    T = eval(Symbol(read(HDF5.attributes(store.file[name])["T"])))
    return Lazy{ProteinStructure{T}}(store.file[name])
end
