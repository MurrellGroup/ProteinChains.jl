function setproperties end
function addproperties end
function removeproperties end

sortnames(np::NamedTuple{names}) where names = NamedTuple{Tuple(sort(collect(names)))}(np)

"""
    AbstractProperty

Abstract type for wrapped properties associated with a [`ProteinChain`](@ref) to define custom behavior.
"""
abstract type AbstractProperty end

const NamedProperties{names} = NamedTuple{names,<:Tuple{Vararg{AbstractProperty}}}

namedproperties(properties::NamedTuple) = map(properties) do value
    value isa AbstractProperty ? value : StandardProperty(value)
end

checkproperty(::Any, ::AbstractProperty) = nothing

unpack(x) = x
unpack(p::AbstractProperty) = p.value

struct StandardProperty{T} <: AbstractProperty
    value::T
end

Base.getindex(p::AbstractProperty, ::Any) = unpack(p)

"""
    IndexableProperty

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
Base.getindex(p::IndexableProperty, i::Union{AbstractVector,Colon}) = selectdim(p.value, ndims(p.value), i) |> IndexableProperty

function checkproperty(parent, p::IndexableProperty)
    if size(p.value, ndims(p.value)) != length(parent)
        throw(DimensionMismatch("Property $(p.value) has length $(size(p.value, ndims(p.value))) but parent has length $(length(parent))"))
    end
    return nothing
end
