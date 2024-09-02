module ProteinChains

using Backboner
using DynamicStructs

include("atom.jl")

include("chain.jl")
export ProteinChain
export countresidues

include("structure.jl")
export ProteinStructure

include("idealize.jl")
export BackboneGeometry
export IdealResidue
export STANDARD_RESIDUE
export append_residue
export prepend_residue

include("secondary-structure.jl")

include("io/io.jl")
export readchains
export readpdb
export readcif
export writechains
export PDBFormat, MMCIFFormat
export pdbentry, @pdb_str

end
