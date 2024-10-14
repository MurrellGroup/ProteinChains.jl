using PeriodicTable: elements
using StaticArrays: SVector

const ELEMENT_SYMBOL_TO_NUMBER = Dict(uppercase(elements[number].symbol) => number for number in 1:118)
const number_TO_ELEMENT_SYMBOL = Dict(n => s for (s, n) in ELEMENT_SYMBOL_TO_NUMBER)

element_symbol_to_number(element_symbol::AbstractString) = ELEMENT_SYMBOL_TO_NUMBER[uppercase(strip(element_symbol))]
number_to_element_symbol(number::Integer) = number_TO_ELEMENT_SYMBOL[number]

function pad_atom_name(name::AbstractString, element_symbol::AbstractString)
    length(name) == 4 && return name
    rpad(" "^(2-length(strip(element_symbol)))*strip(name), 4)
end

encode_atom_name(name::AbstractString, element_symbol::AbstractString) = reinterpret(UInt32, codeunits(pad_atom_name(name, element_symbol)))[1]
decode_atom_name(name::UInt32) = String(reinterpret(UInt8, [name]))

struct Atom{T<:AbstractFloat}
    name::UInt32
    number::Int8
    x::T
    y::T
    z::T
end

Base.convert(::Type{Atom{T}}, atom::Atom) where T = Atom{T}(atom.name, atom.number, atom.x, atom.y, atom.z)

@inline Atom(name::UInt32, number::Integer, x::T, y::T, z::T) where T = Atom{T}(name, number, x, y, z)
@inline Atom(name::AbstractString, element_symbol::AbstractString, x, y, z) =
    Atom(encode_atom_name(name, element_symbol), element_symbol_to_number(element_symbol), x, y, z)

@inline Atom(name, number, coords::AbstractVector) = Atom(name, number, coords...)

atom_name(atom::Atom) = decode_atom_name(atom.name)
atom_number(atom::Atom) = atom.number
atom_coords(atom::Atom) = SVector(atom.x, atom.y, atom.z)
atom_symbol(atom::Atom) = number_to_element_symbol(atom.number)

Base.summary(atom::Atom) = "$(elements[atom.number].name) atom at [$(atom.x), $(atom.y), $(atom.z)])"

Base.show(io::IO, atom::Atom{T}) where T = print(io,
    "Atom(\"$(decode_atom_name(atom.name))\", \"$(number_to_element_symbol(atom.number))\", $(atom_coords(atom)))")