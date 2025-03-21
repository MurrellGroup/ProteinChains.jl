module ProteinChains

using Reexport
@reexport using Backboner
@reexport using AssigningSecondaryStructure
@reexport using DynamicStructs

using Compat: @compat

include("atom.jl")
export Atom
@compat public atom_name, atom_number, atom_coords, atom_symbol, atom_mass

include("properties.jl")
export AbstractProperty, unwrap, Indexable

include("chain.jl")
export ProteinChain
export map_atoms!
export get_atoms, get_backbone
@compat public psi_angles, omega_angles, phi_angles

include("structure.jl")
export ProteinStructure
export map_chains!

include("io/io.jl")
export readcif, readpdb
export writecif, writepdb
export pdbentry, @pdb_str, @mmcifdict_str
export getmmcif, mapmmcif
export BioStructures, MMCIFDict, PDBFormat, MMCIFFormat

include("store.jl")
export ProteinStructureStore
@compat public serialize, deserialize

include("ideal.jl")
export BackboneGeometry, DEFAULT_BACKBONE_GEOMETRY
export IdealResidue, STANDARD_RESIDUE
export append_residue, prepend_residue

include("utils/utils.jl")

include("deprecated/deprecated.jl")

end
