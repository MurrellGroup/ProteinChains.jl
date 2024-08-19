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

Base.propertynames(x::Properties, private::Bool=false) = ((private ? fieldnames(Properties) : ())..., collect(keys(x._dict))...,)

Base.getproperty(x::Properties, name::Symbol) = hasfield(Properties, name) ? getfield(x, name) : getindex(x, name)
Base.setproperty!(x::Properties, name::Symbol, value) = hasfield(Properties, name) ? setfield!(x, name, value) : setindex!(x, value, name)

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

macro properties(expr)
    expr.head != :struct && error("@properties can only be applied to struct definitions")

    struct_mutable = expr.args[1]
    struct_def = expr.args[2]
    struct_body = expr.args[3].args

    if struct_def isa Expr && struct_def.head == :(<:)
        struct_name = struct_def.args[1]
        supertype = struct_def.args[2]
    else
        struct_name = struct_def
        supertype = :Any
    end

    if struct_name isa Expr && struct_name.head == :curly
        type_name = struct_name.args[1]
        type_params = struct_name.args[2:end]
        full_type = Expr(:where, struct_name, type_params...)
    else
        type_name = struct_name
        type_params = []
        full_type = type_name
    end

    push!(struct_body, :(properties::Properties))

    propertynames_expr = quote
        function Base.propertynames(x::$type_name, private::Bool=false)
            return (fieldnames($type_name)..., propertynames(x.properties, false)...)
        end
    end

    getproperty_expr = quote
        function Base.getproperty(x::$type_name, name::Symbol)
            hasfield($type_name, name) && return getfield(x, name)
            name in keys(x.properties) && return getindex(x.properties, name)
            throw(ErrorException("$($type_name) instance has no field or property $name"))
        end
    end

    setproperty_expr = quote
        function Base.setproperty!(x::$type_name, name::Symbol, value)
            hasfield($type_name, name) && return setfield!(x, name, value)
            return setindex!(x.properties, value, name)
        end
    end

    return quote
        $(esc(expr))
        $(esc(propertynames_expr))
        $(esc(getproperty_expr))
        $(esc(setproperty_expr))
    end
end