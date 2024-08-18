using PeriodicTable: elements

const ELEMENT_SYMBOL_TO_ATOMIC_NUMBER = Dict(uppercase(elements[atomic_number].symbol) => atomic_number for atomic_number in 1:118)
const ATOMIC_NUMBER_TO_ELEMENT_SYMBOL = Dict(n => s for (s, n) in ELEMENT_SYMBOL_TO_ATOMIC_NUMBER)

element_symbol_to_atomic_number(element_symbol::AbstractString) = ELEMENT_SYMBOL_TO_ATOMIC_NUMBER[uppercase(element_symbol)]
atomic_number_to_element_symbol(atomic_number::Integer) = ATOMIC_NUMBER_TO_ELEMENT_SYMBOL[atomic_number]

encode_atom_name(atom_name::AbstractString) = reinterpret(UInt32, codeunits(lpad(atom_name, 4)))[1]
decode_atom_name(atom_name::UInt32) = String(reinterpret(UInt8, [atom_name]))

mutable struct Atom{T<:AbstractFloat}
    atom_name::UInt32
    atomic_number::UInt8
    x::T
    y::T
    z::T
end

coords(atom::Atom) = [atom.x, atom.y, atom.z]

Base.summary(atom::Atom) = "$(elements[atom.atomic_number].name) atom at [$(atom.x), $(atom.y), $(atom.z)])"

Atom(atom_name::UInt32, atomic_number::Integer, x::T, y::T, z::T) where T = Atom{T}(atom_name, atomic_number, x, y, z)
Atom(atom_name::UInt32, element_symbol::AbstractString, args...) = Atom(atom_name, element_symbol_to_atomic_number(element_symbol), args...)
Atom(atom_name::AbstractString, args...) = Atom(encode_atom_name(atom_name), args...)
Atom(atom_name::AbstractString, element_symbol, coords::AbstractVector{<:AbstractFloat}) = Atom(atom_name, element_symbol, coords...)
