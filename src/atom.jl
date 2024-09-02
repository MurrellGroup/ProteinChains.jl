using PeriodicTable: elements

const ELEMENT_SYMBOL_TO_ATOMIC_NUMBER = Dict(uppercase(elements[atomic_number].symbol) => atomic_number for atomic_number in 1:118)
const ATOMIC_NUMBER_TO_ELEMENT_SYMBOL = Dict(n => s for (s, n) in ELEMENT_SYMBOL_TO_ATOMIC_NUMBER)

element_symbol_to_atomic_number(element_symbol::AbstractString) = ELEMENT_SYMBOL_TO_ATOMIC_NUMBER[uppercase(strip(element_symbol))]
atomic_number_to_element_symbol(atomic_number::Integer) = ATOMIC_NUMBER_TO_ELEMENT_SYMBOL[atomic_number]

pad_atom_name(atom_name::AbstractString, element_symbol::AbstractString) = rpad(" "^(2-length(strip(element_symbol)))*strip(atom_name), 4)

encode_atom_name(atom_name::AbstractString, element_symbol::AbstractString) = reinterpret(UInt32, codeunits(pad_atom_name(atom_name, element_symbol)))[1]
decode_atom_name(atom_name::UInt32) = String(reinterpret(UInt8, [atom_name]))

mutable struct Atom{T<:AbstractFloat}
    atom_name::UInt32
    atomic_number::Int8
    x::T
    y::T
    z::T
end

Atom(atom_name::UInt32, atomic_number::Integer, x::T, y::T, z::T) where T = Atom{T}(atom_name, atomic_number, x, y, z)
Atom(atom_name::AbstractString, element_symbol::AbstractString, coords::AbstractVector{<:AbstractFloat}) =
    Atom(encode_atom_name(atom_name, element_symbol), element_symbol_to_atomic_number(element_symbol), coords...)

coords(atom::Atom) = [atom.x, atom.y, atom.z]

Base.summary(atom::Atom) = "$(elements[atom.atomic_number].name) atom at [$(atom.x), $(atom.y), $(atom.z)])"

function offset!(atom::Atom, coords::Vector{<:Real})
    @assert length(coords) == 3
    atom.x += coords[1]
    atom.y += coords[2]
    atom.z += coords[3]
    return atom
end 