abstract type AbstractProperty end

Base.getindex(p::AbstractProperty) = p.value

"""
    PersistentProperty(value)

A property of arbitrary type that persists after residue indexing of a chain.

```jldoctest
julia> chain = addproperties(pdb"1ASS"A; x=PersistentProperty(1));

julia> chain.x == chain[1:10].x
true
```
"""
struct PersistentProperty{T} <: AbstractProperty
    value::T
end

Base.getindex(p::PersistentProperty, ::AbstractVector) = p

"""
    IndexableProperty <: AbstractProperty

    IndexableProperty(value::AbstractArray)

An `AbstractArray` property with size `(dims..., length(chain))`, and
residue indexing of the chain being propagated to the last dimension of the array.

```jldoctest
julia> chain = pdb"1ASS"A;

julia> chain = addproperties(pdb"1ASS"A; y=IndexableProperty(rand(2,152)));

julia> chain.y == chain[1:10].y
false

julia> chain.y[:,1:10] == chain[1:10].y
true
```
"""
struct IndexableProperty{T<:AbstractArray} <: AbstractProperty
    value::T
end

Base.getindex(p::IndexableProperty, i::AbstractVector) = selectdim(p.value, ndims(p.value), i) |> IndexableProperty

const NamedProperties{names} = NamedTuple{names,<:Tuple{Vararg{AbstractProperty}}}

sortnames(np::NamedProperties{names}) where names = NamedTuple{Tuple(sort(collect(names)))}(np)

function addproperties end