# # Examples

# ## ProteinStructure and ProteinChain

using ProteinChains

# Fetch a structure from the PDB:

structure = pdb"3NIR"
#-
propertynames(structure)

# We can index by chain ID or by index:

chain = structure["A"]
#-
chain == structure[1]

# A chain has a few basic fields:

propertynames(chain)
#-
chain.id
#-
typeof(chain.atoms)
#-
chain.sequence
#-
chain.numbering
#-
chain.ins_codes # (sometimes just a string of spaces)

# The `atoms` field is a residue-wise vector of vectors of `Atom`s, but often we just want the backbone:

get_backbone(chain)

# Chains can be indexed/reordered with a vector of residue indices:

chain[1:10].numbering
#-
chain[10:-1:1].numbering

# ## Dynamic properties

# The chain type is [*dynamic*](https://github.com/AntonOresten/DynamicStructs.jl), so we can add new properties dynamically at runtime:

chain.taxonomy_id = 3721;
#-
chain.taxonomy_id
#-
delete!(chain, :taxonomy_id);
#-
hasproperty(chain, :taxonomy_id)

# For a property whose last dimension is tied to the number of residues, we can wrap it with `Indexable` to automatically index it when we index the chain:

chain.confidence = Indexable(rand(length(chain)));
#-
chain[1:10].confidence # still wrapped by Indexable
#-
unwrap(chain[1:10].confidence)

# ## Backbone geometry

# There are utility functions for getting the backbone geometry:

get_bond_lengths(chain) # N₁-Ca₁, Ca₁-C₁, C₁-N₂, N₂-Ca₂, ... Caₙ-Cₙ
#-
get_bond_angles(chain) # N₁-Ca₁-C₁, Ca₁-C₁-N₂, C₁-N₂-Ca₂, ... Cₙ-Caₙ-Cₙ
#-
get_torsion_angles(chain) # N₁-Ca₁-C₁-N₂, Ca₁-C₁-N₂-Ca₂, C₁-N₂-Ca₂-C₂, ... Cₙ₋₁-Nₙ-Caₙ-Cₙ
#-

# There are also functions for getting the residue-wise rigid transformation "frames" of the chain:

frames = Frames(chain);
#-
size(frames.rotations)
#-
size(frames.translations)

# ## MMCIF utilities

# MMCIF files contain a lot of information that is not present in PDB files.
# The `mapmmcif` function can be used to map one MMCIF field to another.

mmcifdict = mmcifdict"3HFM";
#-
mapmmcif(mmcifdict, "_atom_site.auth_asym_id"   => "_atom_site.label_entity_id")
#-
mapmmcif(mmcifdict, "_entity_src_gen.entity_id" => "_entity_src_gen.pdbx_gene_src_ncbi_taxonomy_id")
#-
chainid_to_taxonomyid = mapmmcif(mmcifdict,
    "_atom_site.auth_asym_id"   => "_atom_site.label_entity_id",
    "_entity_src_gen.entity_id" => "_entity_src_gen.pdbx_gene_src_ncbi_taxonomy_id")
#-
structure = pdb"3HFM"
for chain in structure
    chain.taxonomy_id = chainid_to_taxonomyid[chain.id]
    println("Set taxonomy_id of chain $(chain.id) to $(chain.taxonomy_id)")
end

# ## ProteinStructureStore

# The `ProteinStructureStore <: AbstractDict{String, ProteinStructure}` type is a lazy wrapper for a file-based storage of `ProteinStructure`s.

dir = mktempdir();
store = ProteinStructureStore(joinpath(dir, "structures.pss"));
#-
store["3NIR"] = pdb"3NIR";
#-
store["3NIR"]
#-
store["3NIR"]["A"]
