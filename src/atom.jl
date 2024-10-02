using PeriodicTable: elements

const ELEMENT_SYMBOL_TO_number = Dict(uppercase(elements[number].symbol) => number for number in 1:118)
const number_TO_ELEMENT_SYMBOL = Dict(n => s for (s, n) in ELEMENT_SYMBOL_TO_number)

element_symbol_to_number(element_symbol::AbstractString) = ELEMENT_SYMBOL_TO_number[uppercase(strip(element_symbol))]
number_to_element_symbol(number::Integer) = number_TO_ELEMENT_SYMBOL[number]

pad_atom_name(name::AbstractString, element_symbol::AbstractString) = rpad(" "^(2-length(strip(element_symbol)))*strip(name), 4)

encode_atom_name(name::AbstractString, element_symbol::AbstractString) = reinterpret(UInt32, codeunits(pad_atom_name(name, element_symbol)))[1]
decode_atom_name(name::UInt32) = String(reinterpret(UInt8, [name]))

struct Atom{T<:AbstractFloat}
    name::UInt32
    number::Int8
    x::T
    y::T
    z::T
end

Atom{T}(atom::Atom) where T = Atom(atom.name, atom.number, T(atom.x), T(atom.y), T(atom.z))

Atom(name::UInt32, number::Integer, x::T, y::T, z::T) where T = Atom{T}(name, number, x, y, z)
Atom(name::AbstractString, element_symbol::AbstractString, coords::AbstractVector{T}) where T =
    Atom(encode_atom_name(name, element_symbol), element_symbol_to_number(element_symbol), coords...)

coords(atom::Atom) = SVector(atom.x, atom.y, atom.z)

Base.summary(atom::Atom) = "$(elements[atom.number].name) atom at [$(atom.x), $(atom.y), $(atom.z)])"

Base.show(io::IO, atom::Atom{T}) where T = print(io,
    "Atom(\"$(decode_atom_name(atom.name))\", \"$(number_to_element_symbol(atom.number))\", $(coords(atom)))")