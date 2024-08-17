module ProteinChains

using Backboner

using BioStructures: BioStructures, PDBFormat, MMCIFFormat

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

include("io.jl")
export readrecord
export readpdb
export readcif

end
