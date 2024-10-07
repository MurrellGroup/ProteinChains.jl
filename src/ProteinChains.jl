module ProteinChains

using StaticArrays: SVector

using Backboner

include("ideal.jl")
export BackboneGeometry
export IdealResidue
export STANDARD_RESIDUE
export append_residue
export prepend_residue

include("atom.jl")
export Atom
export encode_atom_name, decode_atom_name
export coords

include("properties.jl")
export ChainProperty, ResidueProperty

include("chain.jl")
export ProteinChain
export annotate
export map_atoms!
export psi_angles, omega_angles, phi_angles
export get_atoms, get_backbone

include("structure.jl")
export ProteinStructure

include("dataset.jl")
export ProteinDataset

include("io/io.jl")
export readcif, readpdb
export writecif, writepdb
export PDBFormat, MMCIFFormat
export pdbentry, @pdb_str

end
