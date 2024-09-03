module ProteinChains

using DynamicStructs

using Backboner
export Backbone
export ChainedBonds
export get_bond_lengths, get_bond_angles, get_torsion_angles
export Frames

include("ideal.jl")
export BackboneGeometry
export IdealResidue
export STANDARD_RESIDUE
export append_residue
export prepend_residue

include("atom.jl")

include("chain.jl")
export ProteinChain
export countresidues
export psi_angles, omega_angles, phi_angles

include("structure.jl")
export ProteinStructure

include("secondary-structure.jl")

include("io/io.jl")
export readchains
export readpdb
export readcif
export writechains
export PDBFormat, MMCIFFormat
export pdbentry, @pdb_str

end
