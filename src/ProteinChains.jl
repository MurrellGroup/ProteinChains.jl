module ProteinChains

using Backboner

include("properties.jl")

include("chain.jl")
export ProteinChain
export countresidues

include("structure.jl")
export ProteinStructure

include("idealize.jl")
export BackboneGeometry

include("secondary-structure.jl")

include("chain-properties.jl")
export assign_bond_lengths!
export assign_bond_angles!
export assign_torsion_angles!
export assign_secondary_structure!
export assign_standard_residue!
export assign_residue_rotations!
export assign_residue_rotations_quat!
export assign_residue_translations!
export assign_is_knotted!

include("io.jl")
export readchains
export readpdb
export readcif
export writechains
export pdbentry, @pdb_str

end
