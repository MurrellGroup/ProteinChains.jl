abstract type AbstractProperty end

unwrap(x) = x
unwrap(p::AbstractProperty) = p.value

Base.:(==)(p1::T, p2::T) where T<:AbstractProperty = unwrap(p1) == unwrap(p2)

function Base.show(io::IO, ::MIME"text/plain", p::AbstractProperty)
    print(io, typeof(p), " wrapping a ", typeof(p.value).name.name, Base.text_colors[:black], " (get value with `unwrap(x)`)", Base.text_colors[:normal])
end

Base.getindex(p::AbstractProperty, i...) = p

checkproperty(::Any, ::AbstractProperty) = nothing

struct Indexable <: AbstractProperty
    value::AbstractArray
end

Base.getindex(p::Indexable, i) = Indexable(selectdim(p.value, ndims(p.value), i))

checkproperty(parent, p::Indexable) = size(p.value, ndims(p.value)) != length(parent) &&
    throw(DimensionMismatch("Property $(p.value) has length $(size(p.value, ndims(p.value))) but parent has length $(length(parent))"))
