module ProteinChains

using Backboner

using Compat: @compat

include("ideal.jl")
export BackboneGeometry
export IdealResidue, STANDARD_RESIDUE
export append_residue, prepend_residue

include("atom.jl")
export Atom
@compat public (atom_name, atom_number, atom_coords, atom_symbol)

include("properties.jl")
export StandardProperty, IndexableProperty
export addproperties, removeproperties
@compat public AbstractProperty

include("chain.jl")
export ProteinChain
export map_atoms!
export get_atoms, get_backbone
@compat public (psi_angles, omega_angles, phi_angles)

include("structure.jl")
export ProteinStructure
@compat public renumber!

include("io/io.jl")
export readcif, readpdb
export writecif, writepdb
export pdbentry, @pdb_str
export getmmcif, mapmmcif
export BioStructures, MMCIFDict, PDBFormat, MMCIFFormat

include("store/store.jl")
export ProteinStructureStore
@compat public (serialize, deserialize)
@compat public (readattribute, writeattribute)
@compat public (readproperty, writeproperty, deleteproperty)

end
