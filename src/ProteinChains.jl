module ProteinChains

using Backboner

using Compat: @compat

include("atom.jl")
export Atom
@compat public (atom_name, atom_number, atom_coords, atom_symbol, atom_mass)

include("properties.jl")
export AbstractProperty, StandardProperty, IndexableProperty
export setproperties!, addproperties!, removeproperties!
export setproperties, addproperties, removeproperties

include("chain.jl")
export ProteinChain
export map_atoms!
export get_atoms, get_backbone
@compat public (psi_angles, omega_angles, phi_angles)

include("structure.jl")
export ProteinStructure

include("io/io.jl")
export readcif, readpdb
export writecif, writepdb
export pdbentry, @pdb_str, @mmcifdict_str
export getmmcif, mapmmcif
export BioStructures, MMCIFDict, PDBFormat, MMCIFFormat
@compat public renumber

include("store/store.jl")
export ProteinStructureStore
@compat public (serialize, deserialize)
@compat public (readattribute, writeattribute)
@compat public (readproperty, writeproperty, deleteproperty)

end
