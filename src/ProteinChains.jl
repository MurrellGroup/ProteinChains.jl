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
export AbstractProteinChain
export ProteinChain
export length
export psi_angles, omega_angles, phi_angles

include("annotated-chain.jl")
export AnnotatedProteinChain
export annotate!, annotate
export annotate_indexable!, annotate_indexable

include("structure.jl")
export ProteinStructure

include("secondary-structure.jl")
export assign_secondary_structure, assign_secondary_structure!

include("io/io.jl")
export readcif, readpdb
export writecif, writepdb
export PDBFormat, MMCIFFormat
export pdbentry, @pdb_str

end
