using PeriodicTable: elements
using StaticArrays: SVector

const ELEMENT_SYMBOL_TO_NUMBER = Dict(uppercase(elements[number].symbol) => number for number in 1:118)
const NUMBER_TO_ELEMENT_SYMBOL = Dict(n => s for (s, n) in ELEMENT_SYMBOL_TO_NUMBER)

element_symbol_to_number(element_symbol::AbstractString) = get(ELEMENT_SYMBOL_TO_NUMBER, uppercase(strip(element_symbol)), 0)
number_to_element_symbol(number::Integer) = get(NUMBER_TO_ELEMENT_SYMBOL, number, "X")

function pad_atom_name(name::AbstractString, element_symbol::AbstractString)
    length(name) == 4 && return name
    rpad(" "^(2-length(strip(element_symbol)))*strip(name), 4)
end

using BitIntegers

BitIntegers.@define_integers 24

primitive type AtomName 24 end

function Base.convert(::Type{AtomName}, name::AbstractString)
    x = zero(UInt24)
    for (i, c) in enumerate(codeunits(name))
        x |= ((c - 0x20) % UInt24) << (6 * (i - 1))
    end
    reinterpret(AtomName, x)
end

function Base.convert(::Type{String}, name::AtomName)
    x = reinterpret(UInt24, name)
    cs = Vector{UInt8}(undef, 4)
    for i in 1:4
        cs[i] = x & 0b111111 + 0x20
        x >>= 6
    end
    String(cs)
end

encode_atom_name(name::AbstractString, element_symbol::AbstractString) = convert(AtomName, pad_atom_name(name, element_symbol))
decode_atom_name(name::AtomName) = convert(String, name)

struct Atom{T<:Real}
    name::AtomName
    number::Int8
    x::T
    y::T
    z::T
end

Base.convert(::Type{Atom{T}}, atom::Atom) where T = Atom{T}(atom.name, atom.number, atom.x, atom.y, atom.z)

@inline Atom(name::AtomName, number::Integer, x::T, y::T, z::T) where T = Atom{T}(name, number, x, y, z)
@inline Atom(name::AbstractString, element_symbol::AbstractString, x, y, z) =
    Atom(encode_atom_name(name, element_symbol), element_symbol_to_number(element_symbol), x, y, z)

@inline Atom(name, number, coords::AbstractVector) = Atom(name, number, coords...)

atom_name(atom::Atom) = decode_atom_name(atom.name)
atom_number(atom::Atom) = atom.number
atom_coords(atom::Atom) = SVector(atom.x, atom.y, atom.z)
atom_symbol(atom::Atom) = number_to_element_symbol(atom.number)
atom_mass(atom::Atom) = elements[atom_number(atom)].atomic_mass

Base.summary(atom::Atom) = "$(elements[atom.number].name) atom at [$(atom.x), $(atom.y), $(atom.z)])"

Base.show(io::IO, atom::Atom{T}) where T = print(io,
    "Atom(\"$(decode_atom_name(atom.name))\", \"$(number_to_element_symbol(atom.number))\", $(atom_coords(atom)))")
