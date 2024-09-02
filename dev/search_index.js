var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = ProteinChains","category":"page"},{"location":"#ProteinChains","page":"Home","title":"ProteinChains","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for ProteinChains.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [ProteinChains]","category":"page"},{"location":"#ProteinChains.STANDARD_RESIDUE_TEMPLATE","page":"Home","title":"ProteinChains.STANDARD_RESIDUE_TEMPLATE","text":"STANDARD_RESIDUE_TEMPLATE\n\nThis is a template of a \"standard residue\", with a very specific and distinct shape, size, and orientation. which needs to be consistent if we want to represent protein structures as sets of residue rotations and translations.\n\nThus, we can use this residue as a template for aligning other residues with very precise geometry to it.\n\njulia> IdealResidue{Float64}(BackboneGeometry(N_Ca_C_angle = 1.93); template=ProteinChains.STANDARD_RESIDUE)\n3×3 IdealResidue{Float64}:\n -1.06447   -0.199174   1.26364\n  0.646303  -0.529648  -0.116655\n  0.0        0.0        0.0\n\n\n\n\n\n","category":"constant"},{"location":"#ProteinChains.BackboneGeometry","page":"Home","title":"ProteinChains.BackboneGeometry","text":"BackboneGeometry(;\n    N_Ca_length = 1.46,\n    Ca_C_length = 1.52,\n    C_N_length = 1.33,\n\n    N_Ca_C_angle = 1.94,\n    Ca_C_N_angle = 2.03,\n    C_N_Ca_angle = 2.13,\n)\n\nDefine the idealized bond lengths and bond angles of a protein backbone.\n\n\n\n\n\n","category":"type"},{"location":"#ProteinChains.IdealResidue","page":"Home","title":"ProteinChains.IdealResidue","text":"IdealResidue{T<:AbstractFloat} <: AbstractMatrix{T}\n\nIdealResidue{T}(backbone_geometry=DEFAULT_BACKBONE_GEOMETRY; template=nothing) where T\n\nA 3x3 matrix representing the idealized geometry of a protein residue, with columns representing the N, Ca, and C atom positions of a residue positioned at the origin.\n\n\n\n\n\n","category":"type"},{"location":"#ProteinChains.ProteinChain","page":"Home","title":"ProteinChains.ProteinChain","text":"ProteinChain{T<:AbstractFloat}\n\nExamples\n\njulia> structure = pdb\"1EYE\";\n[ Info: Downloading file from PDB: 1EYE\n\njulia> structure[1]\n253-residue ProteinChain \"A\":\n  4 fields:\n    id::String = \"A\"\n    sequence::String = <field value exceeds max length>\n    backbone::Array{Float64,3} = <field value exceeds max length>\n    atoms::Vector{Vector{ProteinChains.Atom{Float64}}} = <field value exceeds max length>\n  2 properties:\n    numbering::Vector{Int64} = [5, 6, 7, 8, 9, 10, 11, 12, 13, 14  …  265, 266, 267, 268, 269, 270, 271, 272, 273, 274]\n    modelnum::Int64 = 1\n\n\n\n\n\n","category":"type"},{"location":"#ProteinChains.ProteinStructure","page":"Home","title":"ProteinChains.ProteinStructure","text":"ProteinStructure{T<:AbstractFloat} <: AbstractVector{ProteinChain{T}}\n\nExamples\n\njulia> structure = pdb\"1EYE\"\n1-chain ProteinStructure \"1EYE.cif\" with 2 dynamic properties:\n  2 fields:\n    name::String = \"1EYE.cif\"\n    chains::Vector{ProteinChain{Float64}} = <field value exceeds max length>\n  2 properties:\n    ids::Vector{String} = [\"A\"]\n    lengths::Vector{Int64} = [253]\n\njulia> structure[\"A\"] isa ProteinChain\ntrue\n\n\n\n\n\n","category":"type"},{"location":"#ProteinChains.append_residue-Tuple{Backboner.Backbone, Vector{<:Real}}","page":"Home","title":"ProteinChains.append_residue","text":"append_residue(Backbone::Backbone, torsion_angles::Vector{<:Real}; ideal::BackboneGeometry=DEFAULT_BACKBONE_GEOMETRY)\n\nCreate a new backbone by appending 3 new torsion angles (ψ, ω, ϕ) at the end, using bond lengths and bond angles specified in BackboneGeometry.\n\n\n\n\n\n","category":"method"},{"location":"#ProteinChains.prepend_residue-Tuple{Backboner.Backbone, Vector{<:Real}}","page":"Home","title":"ProteinChains.prepend_residue","text":"append_residue(Backbone::Backbone, torsion_angles::Vector{<:Real}; ideal::BackboneGeometry=DEFAULT_BACKBONE_GEOMETRY)\n\nCreate a new backbone by prepending 3 new torsion angles (ψ, ω, ϕ) at the beginning, using bond lengths and bond angles specified in the BackboneGeometry.\n\nnote: Note\nThe torsion angle order is the same as it would be when appending. The order is not reversed.\n\n\n\n\n\n","category":"method"},{"location":"#ProteinChains.renumber!-Tuple{ProteinStructure, BioStructures.MMCIFDict}","page":"Home","title":"ProteinChains.renumber!","text":"renumber!(structure::ProteinStructure, mmcif_dict::BioStructures.MMCIFDict)\n\nRenumber the residues in a ProteinStructure object according to the numbering aligned to a reference sequence in the MMCIF file.\n\n\n\n\n\n","category":"method"}]
}