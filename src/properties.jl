function _getproperty(pc::T, property::Symbol) where T
    property in fieldnames(T) && return getfield(pc, property)
    property in keys(pc.properties) && return pc.properties[property]
    throw(ErrorException("$(T) instance has no field or property $property"))
end

function _setproperty!(pc::T, property::Symbol, x) where T
    property in fieldnames(T) && return setfield!(pc, property, x)
    return pc.properties[property] = x
end

truncate(s::AbstractString, len::Integer) = length(s) > len ? String(first(s, len-1) * '…') : s

function printfield(io::IO, x, field::Symbol; indent=4)
    value = getproperty(x, field)
    print(io, "\n"*" "^indent)
    printstyled(io, string(field); color=:white)
    printstyled(io, "::"; color=:red)
    printstyled(io, string(typeof(value)); color=:blue)
    printstyled(io, " = "; color=:red)
    printstyled(io, truncate(repr(value; context=io), 100))
end

function printproperty(io::IO, x, property::Symbol; indent=4)
    value = getproperty(x, property)
    print(io, "\n"*" "^indent)
    printstyled(io, ":", string(property); color=:magenta)
    printstyled(io, " => "; color=:red)
    printstyled(io, truncate(repr(value; context=io), 100))
end

# TODO: `assign_property!(chain, property)` vs `assign_<property>!(chain)`
#assignproperty!(x, property::Symbol, args...) = assignproperty!(x, Val(property), args...)