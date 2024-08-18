using OrderedCollections: OrderedDict

export Properties

struct Properties <: AbstractDict{Symbol,Any}
    _dict::OrderedDict{Symbol,Any}
end

Properties() = Properties(OrderedDict{Symbol,Any}())
Properties(pairs::Vararg{<:Pair{Symbol}}) = Properties(OrderedDict(name => value for (name, value) in pairs))
Properties(args...; kwargs...) = Properties(args..., kwargs...)

Base.length(x::Properties) = length(x._dict)
Base.iterate(x::Properties, args...) = iterate(x._dict, args...)
Base.keys(x::Properties) = keys(x._dict)
Base.values(x::Properties) = values(x._dict)

Base.getindex(x::Properties, name::Symbol) = getindex(x._dict, name)
Base.setindex!(x::Properties, value, name::Symbol) = setindex!(x._dict, value, name)

Base.hasproperty(x::Properties, name::Symbol) = hasfield(Properties, name) || name in keys(x)
Base.getproperty(x::Properties, name::Symbol) = hasfield(Properties, name) ? getfield(x, name) : getindex(x, name)
Base.setproperty!(x::Properties, name::Symbol, value) = hasfield(Properties, name) ? setfield!(x, name, value) : setindex!(x, value, name)

Base.propertynames(x::Properties, private::Bool=false) = ((private ? fieldnames(Properties) : ())..., collect(keys(x._dict))...,)

struct HasProperties{T}
    _parent::T

    function HasProperties(_parent::T) where T
        hasfield(T, :properties) && getfield(_parent, :properties) isa Properties || throw(ArgumentError("input type `$T` needs to have a properties::$(Properties) field"))
        new{T}(_parent)
    end
end

Base.hasproperty(x::HasProperties{T}, name::Symbol) where T = hasfield(T, name) || hasfield(HasProperties, name) || name in keys(x._parent.properties)

function Base.getproperty(x::HasProperties{T}, name::Symbol) where T
    hasfield(HasProperties, name) && return getfield(x, name)
    hasfield(T, name) && return getfield(x._parent, name)
    name in keys(x._parent.properties) && return getindex(x._parent.properties, name)
    throw(ErrorException("$(T) instance has no field or property $name"))
end

function Base.setproperty!(x::HasProperties{T}, name::Symbol, value) where T
    hasfield(HasProperties, name) && return setfield(x, name, value)
    hasfield(T, name) && return setfield!(x._parent, name, value)
    setindex!(x._parent.properties, value, name)
end

Base.propertynames(x::HasProperties{T}, private::Bool=false) where T = ((private ? fieldnames(HasProperties) : ())..., fieldnames(T)..., propertynames(x._parent.properties, false)...)

truncate(s::AbstractString, len::Integer) = length(s) > len ? String(first(s, len-1) * '…') : s

function printfield(io::IO, x, name::Symbol; indent=4)
    value = getproperty(x, name)
    print(io, "\n"*" "^indent)
    printstyled(io, string(name); color=:white)
    printstyled(io, "::"; color=:red)
    printstyled(io, replace(string(typeof(value)), " " => ""); color=:blue)
    printstyled(io, " = "; color=:red)
    printstyled(io, truncate(repr(value; context=io), 100))
end

function showproperties(io::IO, x::T) where T
    context = IOContext(io, :compact => true, :limit => true)
    print(context, summary(x), ":")
    printstyled(context, "\n  fields:", color=:yellow)
    for fieldname in fieldnames(T)[1:end-1]
        printfield(context, x, fieldname)
    end
    printstyled(context, "\n  properties:", color=:yellow)
    isempty(x.properties) && print(io, " (none)")
    for name in keys(x.properties)
        printfield(context, x, name)
    end
end
