abstract type AbstractProperty end

Base.getindex(p::AbstractProperty) = p.value

struct ChainProperty{T} <: AbstractProperty
    value::T
end

Base.getindex(prop::ChainProperty, ::AbstractVector) = prop

struct ResidueProperty{T<:AbstractArray} <: AbstractProperty
    value::T
end

Base.getindex(prop::ResidueProperty, i::AbstractVector) = selectdim(prop.value, ndims(prop.value), i) |> ResidueProperty

const NamedProperties{names} = NamedTuple{names,<:Tuple{Vararg{AbstractProperty}}}

function annotate end