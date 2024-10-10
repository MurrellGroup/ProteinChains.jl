module ProteinChains

using Backboner

using Compat: @compat

include("ideal.jl")
export BackboneGeometry
export IdealResidue, STANDARD_RESIDUE
export append_residue, prepend_residue

include("atom.jl")
export Atom
@compat public (atom_name, atom_number, atom_coords)

include("properties.jl")
export PersistentProperty, IndexableProperty
export addproperties
@compat public (AbstractProperty, NamedProperties)

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
export pdbentry, @pdb_str
export PDBFormat, MMCIFFormat

include("store/store.jl")
export ProteinStructureStore
@compat public (serialize, deserialize)

end
