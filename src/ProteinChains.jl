module ProteinChains

using Backboner
using DynamicStructs

using Compat: @compat

include("atom.jl")
export Atom
@compat public (atom_name, atom_number, atom_coords, atom_symbol, atom_mass)

include("properties.jl")
export AbstractProperty, unwrap, Indexable

include("chain.jl")
export ProteinChain
export map_atoms!
export get_atoms, get_backbone
@compat public (psi_angles, omega_angles, phi_angles)

include("structure.jl")
export ProteinStructure
export map_chains!

include("io/io.jl")
export readcif, readpdb
export writecif, writepdb
export pdbentry, @pdb_str, @mmcifdict_str
export getmmcif, mapmmcif
export BioStructures, MMCIFDict, PDBFormat, MMCIFFormat

include("store/store.jl")
export ProteinStructureStore
@compat public (serialize, deserialize)
@compat public (readattribute, writeattribute)
@compat public (readproperty, writeproperty, deleteproperty)

include("ideal.jl")
export BackboneGeometry, DEFAULT_BACKBONE_GEOMETRY
export IdealResidue, STANDARD_RESIDUE
export append_residue, prepend_residue

using PrecompileTools

@compile_workload begin
    read("src/precompile/3NIR.pdb", ProteinStructure)
    structure = read("src/precompile/3NIR.cif", ProteinStructure)
    io = IOBuffer()
    show(io, MIME("text/plain"), structure)
end

end
