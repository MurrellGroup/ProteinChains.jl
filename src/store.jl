import JLD2
using InlineStrings: String31

"""
    ProteinStructureStore <: AbstractDict{InlineStrings.String31,ProteinStructure}

A JLD2-based store for protein structures implementing the `AbstractDict` interface,
allowing for dictionary operations on the stored structures.

Keys are stored as `InlineStrings.String31` objects to reduce references.
This means keys are limited to 31 bytes.

A `ProteinStructureStore` gets closed automatically when there no longer exists a program-accessible reference to it.

## Examples

```jldoctest
julia> store = ProteinStructureStore("store.pss")
ProteinStructureStore with 0 entries

julia> store["3HFM"] = pdb"3HFM"
[ Info: Downloading file from PDB: 3HFM
3-element ProteinStructure "3HFM.cif":
 215-residue ProteinChain{Float64} (H)
 214-residue ProteinChain{Float64} (L)
 129-residue ProteinChain{Float64} (Y)

julia> store
ProteinStructureStore with 1 entry

julia> keys(store)
Set{InlineStrings.String31} with 1 element:
  InlineStrings.String31("3HFM")

julia> delete!(store, "3HFM")
ProteinStructureStore with 0 entries
```
"""
mutable struct ProteinStructureStore <: AbstractDict{String,ProteinStructure}
    filename::String
    file::JLD2.JLDFile
    keys::Set{String31}
    mode::String
end

function ProteinStructureStore(filename::AbstractString, mode::AbstractString="a+")
    file = JLD2.jldopen(filename, mode)
    store = ProteinStructureStore(filename, file, Set(keys(file)), mode)
    finalizer(close, store)
    return store
end


Base.close(store::ProteinStructureStore) = close(store.file)

Base.keys(store::ProteinStructureStore) = store.keys
Base.length(store::ProteinStructureStore) = length(store.keys)
Base.getindex(store::ProteinStructureStore, key::AbstractString) = store.file[key]
Base.get(store::ProteinStructureStore, key, default) = key in keys(store) ? store[key] : default


function Base.delete!(store::ProteinStructureStore, key::AbstractString)
    if haskey(store, key)
        delete!(store.file, key)
        delete!(store.keys, key)
    end
    return store
end


function Base.setindex!(store::ProteinStructureStore, value::ProteinStructure, key::AbstractString)
    haskey(store, key) && delete!(store, key)
    store.file[key] = value
    push!(store.keys, key)
    return store
end

function Base.iterate(store::ProteinStructureStore, state...)
    iterate_step((key,state)::Tuple{AbstractString,Int}) = (key => store[key], state)
    iterate_step(::Nothing) = nothing
    return iterate_step(iterate(keys(store), state...))
end

Base.show(io::IO, store::ProteinStructureStore) = print(io, "$ProteinStructureStore(\"$(store.filename)\", $(store.mode))")
Base.show(io::IO, ::MIME"text/plain", store::ProteinStructureStore) = print(io, summary(store))

Base.open(::Type{ProteinStructureStore}, filename::AbstractString, args...) = ProteinStructureStore(filename, args...)

"""
    ProteinStructureStore(f::Function, filename, mode="a+")
"""
ProteinStructureStore(f::Function, args...) = open(f, ProteinStructureStore, args...)

"""
    serialize(filename::AbstractString, structures::AbstractVector{<:ProteinStructure})

Serialize a vector of `ProteinStructure` objects to a JLD2 file.
This function creates a new `ProteinStructureStore` and writes each structure in the input vector to it.
Each structure is stored using its name as the key.
"""
function serialize(filename::AbstractString, structures::AbstractVector{<:ProteinStructure})
    ProteinStructureStore(filename, "a+") do store
        for structure in structures
            store[structure.name] = structure
        end
    end
    return filename
end

"""
    deserialize(filename::AbstractString)

Deserialize `ProteinStructure` objects from a JLD2 file.
Returns a `Vector{ProteinStructure}` of all structures stored in the file.
"""
deserialize(filename::AbstractString) = ProteinStructureStore(collect ∘ values, filename, "r")
