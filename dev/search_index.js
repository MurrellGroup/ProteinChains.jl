var documenterSearchIndex = {"docs":
[{"location":"generated/examples/","page":"Examples","title":"Examples","text":"EditURL = \"../../../examples/examples.jl\"","category":"page"},{"location":"generated/examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"generated/examples/#ProteinStructure-and-ProteinChain","page":"Examples","title":"ProteinStructure and ProteinChain","text":"","category":"section"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"using ProteinChains","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"Fetch a structure from the PDB:","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"structure = pdb\"3NIR\"","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"propertynames(structure)","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"We can index by chain ID or by index:","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"chain = structure[\"A\"]","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"chain == structure[1]","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"A chain has a few basic fields:","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"propertynames(chain)","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"chain.id","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"typeof(chain.atoms)","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"chain.sequence","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"chain.numbering","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"chain.ins_codes # (sometimes just a string of spaces)","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"The atoms field is a residue-wise vector of vectors of Atoms, but often we just want the backbone:","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"get_backbone(chain)","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"Chains can be indexed/reordered with a vector of residue indices:","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"chain[1:10].numbering","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"chain[10:-1:1].numbering","category":"page"},{"location":"generated/examples/#Dynamic-properties","page":"Examples","title":"Dynamic properties","text":"","category":"section"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"The chain type is dynamic, so we can add new properties dynamically at runtime:","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"chain.taxonomy_id = 3721;\nnothing #hide","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"chain.taxonomy_id","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"delete!(chain, :taxonomy_id);\nnothing #hide","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"hasproperty(chain, :taxonomy_id)","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"For a property whose last dimension is tied to the number of residues, we can wrap it with Indexable to automatically index it when we index the chain:","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"chain.confidence = Indexable(rand(length(chain)));\nnothing #hide","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"chain[1:10].confidence # still wrapped by Indexable","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"unwrap(chain[1:10].confidence)","category":"page"},{"location":"generated/examples/#Backbone-geometry","page":"Examples","title":"Backbone geometry","text":"","category":"section"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"There are utility functions for getting the backbone geometry:","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"get_bond_lengths(chain) # N₁-Ca₁, Ca₁-C₁, C₁-N₂, N₂-Ca₂, ... Caₙ-Cₙ","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"get_bond_angles(chain) # N₁-Ca₁-C₁, Ca₁-C₁-N₂, C₁-N₂-Ca₂, ... Cₙ-Caₙ-Cₙ","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"get_torsion_angles(chain) # N₁-Ca₁-C₁-N₂, Ca₁-C₁-N₂-Ca₂, C₁-N₂-Ca₂-C₂, ... Cₙ₋₁-Nₙ-Caₙ-Cₙ","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"There are also functions for getting the residue-wise rigid transformation \"frames\" of the chain:","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"frames = Frames(chain);\nnothing #hide","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"size(frames.rotations)","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"size(frames.translations)","category":"page"},{"location":"generated/examples/#MMCIF-utilities","page":"Examples","title":"MMCIF utilities","text":"","category":"section"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"MMCIF files contain a lot of information that is not present in PDB files. The mapmmcif function can be used to map one MMCIF field to another.","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"mmcifdict = mmcifdict\"3HFM\";\nnothing #hide","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"mapmmcif(mmcifdict, \"_atom_site.auth_asym_id\"   => \"_atom_site.label_entity_id\")","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"mapmmcif(mmcifdict, \"_entity_src_gen.entity_id\" => \"_entity_src_gen.pdbx_gene_src_ncbi_taxonomy_id\")","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"chainid_to_taxonomyid = mapmmcif(mmcifdict,\n    \"_atom_site.auth_asym_id\"   => \"_atom_site.label_entity_id\",\n    \"_entity_src_gen.entity_id\" => \"_entity_src_gen.pdbx_gene_src_ncbi_taxonomy_id\")","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"structure = pdb\"3HFM\"\nfor chain in structure\n    chain.taxonomy_id = chainid_to_taxonomyid[chain.id]\n    println(\"Set taxonomy_id of chain $(chain.id) to $(chain.taxonomy_id)\")\nend","category":"page"},{"location":"generated/examples/#ProteinStructureStore","page":"Examples","title":"ProteinStructureStore","text":"","category":"section"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"The ProteinStructureStore <: AbstractDict{String, ProteinStructure} type is a lazy wrapper for a file-based storage of ProteinStructures.","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"dir = mktempdir();\nstore = ProteinStructureStore(joinpath(dir, \"structures.pss\"));\nnothing #hide","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"store[\"3NIR\"] = pdb\"3NIR\";\nnothing #hide","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"store[\"3NIR\"]","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"store[\"3NIR\"][\"A\"]","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"","category":"page"},{"location":"generated/examples/","page":"Examples","title":"Examples","text":"This page was generated using Literate.jl.","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = ProteinChains","category":"page"},{"location":"#ProteinChains","page":"Home","title":"ProteinChains","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for ProteinChains.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [ProteinChains]","category":"page"},{"location":"#ProteinChains.STANDARD_RESIDUE_TEMPLATE","page":"Home","title":"ProteinChains.STANDARD_RESIDUE_TEMPLATE","text":"STANDARD_RESIDUE_TEMPLATE\n\nThis is a template of a \"standard residue\", with a very specific and distinct shape, size, and orientation. which needs to be consistent if we want to represent protein structures as sets of residue rotations and translations.\n\nThus, we can use this residue as a template for aligning other residues with very precise geometry to it.\n\njulia> IdealResidue{Float64}(BackboneGeometry(N_Ca_C_angle = 1.93); template=ProteinChains.STANDARD_RESIDUE_TEMPLATE)\n3×3 IdealResidue{Float64}:\n -1.06447   -0.199174   1.26364\n  0.646303  -0.529648  -0.116655\n  0.0        0.0        0.0\n\n\n\n\n\n","category":"constant"},{"location":"#ProteinChains.BackboneGeometry","page":"Home","title":"ProteinChains.BackboneGeometry","text":"BackboneGeometry(;\n    N_Ca_length = 1.46,\n    Ca_C_length = 1.52,\n    C_N_length = 1.33,\n\n    N_Ca_C_angle = 1.94,\n    Ca_C_N_angle = 2.03,\n    C_N_Ca_angle = 2.13,\n)\n\nDefine the idealized bond lengths and bond angles of a protein backbone.\n\n\n\n\n\n","category":"type"},{"location":"#ProteinChains.IdealResidue","page":"Home","title":"ProteinChains.IdealResidue","text":"IdealResidue{T<:AbstractFloat} <: AbstractMatrix{T}\n\nIdealResidue{T}(backbone_geometry=DEFAULT_BACKBONE_GEOMETRY; template=nothing) where T\n\nA 3x3 matrix representing the idealized geometry of a protein residue, with columns representing the N, Ca, and C atom positions of a residue positioned at the origin.\n\n\n\n\n\n","category":"type"},{"location":"#ProteinChains.ProteinChain","page":"Home","title":"ProteinChains.ProteinChain","text":"ProteinChain{T<:Real}\n\nRepresents a protein chain with a basic set of fields from which some other properties might be derived.\n\n\n\n\n\n","category":"type"},{"location":"#ProteinChains.ProteinStructure","page":"Home","title":"ProteinChains.ProteinStructure","text":"ProteinStructure{T} <: AbstractVector{ProteinChain{T}}\n\nFields\n\nname::String: Usually just the base name of the original file.\natoms::Vector{Atom{T}}: free atoms from the structure that were not part of any protein residue.\nchains::Vector{ProteinChain{T}}: a collection of ProteinChains.\n\n\n\n\n\n","category":"type"},{"location":"#ProteinChains.ProteinStructure-Union{Tuple{T}, Tuple{AbstractString, Array{ProteinChain{T}, 1}}} where T","page":"Home","title":"ProteinChains.ProteinStructure","text":"ProteinStructure(name, chains; properties...)\n\n\n\n\n\n","category":"method"},{"location":"#ProteinChains.ProteinStructureStore","page":"Home","title":"ProteinChains.ProteinStructureStore","text":"ProteinStructureStore <: AbstractDict{String,ProteinStructure}\n\nAn JLD2-based store for protein structures implementing the AbstractDict interface, allowing for dictionary operations on the stored structures.\n\nA ProteinStructureStore gets closed automatically when there no longer exists a program-accessible reference to it.\n\nExamples\n\njulia> store = ProteinStructureStore(\"store.pss\")\nProteinStructureStore with 0 entries\n\njulia> store[\"3HFM\"] = pdb\"3HFM\"\n[ Info: Downloading file from PDB: 3HFM\n3-element ProteinStructure \"3HFM.cif\":\n 215-residue ProteinChain{Float64} (H)\n 214-residue ProteinChain{Float64} (L)\n 129-residue ProteinChain{Float64} (Y)\n\njulia> store\nProteinStructureStore with 1 entry\n\njulia> keys(store)\nSet{String} with 1 element:\n  \"3HFM\"\n\njulia> delete!(store, \"3HFM\")\nProteinStructureStore with 0 entries\n\n\n\n\n\n","category":"type"},{"location":"#ProteinChains.ProteinStructureStore-Tuple{Function, Vararg{Any}}","page":"Home","title":"ProteinChains.ProteinStructureStore","text":"ProteinStructureStore(f::Function, filename, mode=\"a+\")\n\n\n\n\n\n","category":"method"},{"location":"#ProteinChains.append_residue-Tuple{Backbone, Vector{<:Real}}","page":"Home","title":"ProteinChains.append_residue","text":"append_residue(Backbone::Backbone, torsion_angles::Vector{<:Real}; ideal::BackboneGeometry=DEFAULT_BACKBONE_GEOMETRY)\n\nCreate a new backbone by appending 3 new torsion angles (ψ, ω, ϕ) at the end, using bond lengths and bond angles specified in BackboneGeometry.\n\n\n\n\n\n","category":"method"},{"location":"#ProteinChains.deserialize-Tuple{AbstractString}","page":"Home","title":"ProteinChains.deserialize","text":"deserialize(filename::AbstractString)\n\nDeserialize ProteinStructure objects from a JLD2 file. Returns a Vector{ProteinStructure} of all structures stored in the file.\n\n\n\n\n\n","category":"method"},{"location":"#ProteinChains.mapmmcif-Tuple{Any, Vararg{Pair{String, String}}}","page":"Home","title":"ProteinChains.mapmmcif","text":"mapmmcif(mmcifdict, field1 => field2, field3 => field4, ...)\n\njulia> import BioStructures\n\njulia> filename = BioStructures.downloadpdb(\"3HFM\", format=BioStructures.MMCIFFormat);\n[ Info: Downloading file from PDB: 3HFM\n\njulia> mmcifdict = BioStructures.MMCIFDict(filename);\n\njulia> mapmmcif(mmcifdict,\n           \"_atom_site.auth_asym_id\"   => \"_atom_site.label_entity_id\",\n           \"_entity_src_gen.entity_id\" => \"_entity_src_gen.pdbx_gene_src_ncbi_taxonomy_id\")\nDict{String, String} with 3 entries:\n  \"Y\" => \"9031\"\n  \"L\" => \"10090\"\n  \"H\" => \"10090\"\n\n\n\n\n\n","category":"method"},{"location":"#ProteinChains.pdbentry-Tuple{AbstractString}","page":"Home","title":"ProteinChains.pdbentry","text":"pdbentry(pdbid::AbstractString; format=MMCIFFormat, kwargs...)\n\nKeyword arguments get propagated to BioStructures.downloadpdb\n\nDownloads are cached in a temporary directory.\n\nExamples\n\njulia> pdbentry(\"1EYE\")\n[ Info: Downloading file from PDB: 1EYE\n1-element ProteinStructure \"1EYE.cif\":\n 256-residue ProteinChain{Float64} (A)\n\njulia> pdb\"1EYE\" # string macro for convenience\n[ Info: File exists: 1EYE\n1-element ProteinStructure \"1EYE.cif\":\n 256-residue ProteinChain{Float64} (A)\n\njulia> pdb\"1EYE\"A # string suffix to get a specific chain\n[ Info: File exists: 1EYE\n256-residue ProteinChain{Float64} (A)\n\njulia> pdb\"1EYE\"1 # integer suffix to specify \"ba_number\" keyword\n[ Info: Downloading file from PDB: 1EYE\n2-element ProteinStructure \"1EYE_ba1.cif\":\n 256-residue ProteinChain{Float64} (A)\n 256-residue ProteinChain{Float64} (A-2)\n\n\n\n\n\n","category":"method"},{"location":"#ProteinChains.prepend_residue-Tuple{Backbone, Vector{<:Real}}","page":"Home","title":"ProteinChains.prepend_residue","text":"append_residue(Backbone::Backbone, torsion_angles::Vector{<:Real}; ideal::BackboneGeometry=DEFAULT_BACKBONE_GEOMETRY)\n\nCreate a new backbone by prepending 3 new torsion angles (ψ, ω, ϕ) at the beginning, using bond lengths and bond angles specified in the BackboneGeometry.\n\nnote: Note\nThe torsion angle order is the same as it would be when appending. The order is not reversed.\n\n\n\n\n\n","category":"method"},{"location":"#ProteinChains.serialize-Tuple{AbstractString, AbstractVector{<:ProteinStructure}}","page":"Home","title":"ProteinChains.serialize","text":"serialize(filename::AbstractString, structures::AbstractVector{<:ProteinStructure})\n\nSerialize a vector of ProteinStructure objects to a JLD2 file. This function creates a new ProteinStructureStore and writes each structure in the input vector to it. Each structure is stored using its name as the key.\n\n\n\n\n\n","category":"method"}]
}
