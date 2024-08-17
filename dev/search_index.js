var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = ProteinChains","category":"page"},{"location":"#ProteinChains","page":"Home","title":"ProteinChains","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for ProteinChains.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [ProteinChains]","category":"page"},{"location":"#ProteinChains.ProteinChain","page":"Home","title":"ProteinChains.ProteinChain","text":"ProteinChain\n\nExamples\n\njulia> structure = pdb\"1EYE\"\n[ Info: Downloading file from PDB: 1EYE\n1-chain ProteinStructure \"1EYE.cif\" with 0 properties:\n  256-residue ProteinChain \"A\" with 0 properties\n\njulia> structure[1]\n256-residue ProteinChain \"A\" with 0 properties:\n  fields:\n    id::String = \"A\"\n    aminoacids::String = \"PVQVMGVLNVTDDSFSDGGCYLDLDDAVKHGLAMAAAGAGIVDVGGETSRVIPVVKELAAQGITVSIDTMRADVARAALQNGAQMVNDVSGGRADPAM…\n    backbone::Array{Float64, 3} = [45.592 44.171 43.719; -10.864 -10.936 -9.688; 30.192 30.504 31.278;;; 42.568 42.02 40.707; -9.163 …\n    numbers::Vector{Int64} = [5, 6, 7, 8, 9, 10, 11, 12, 13, 14  …  265, 266, 267, 268, 269, 270, 271, 272, 273, 274]\n  properties: (none)\n\n\n\n\n\n","category":"type"},{"location":"#ProteinChains.ProteinStructure","page":"Home","title":"ProteinChains.ProteinStructure","text":"ProteinStructure <: AbstractVector{ProteinChain}\n\n\n\n\n\n","category":"type"},{"location":"#ProteinChains.readcif-Tuple{AbstractString}","page":"Home","title":"ProteinChains.readcif","text":"readcif(path) -> chains::Vector{ProteinChain}\n\n\n\n\n\n","category":"method"},{"location":"#ProteinChains.readpdb-Tuple{AbstractString}","page":"Home","title":"ProteinChains.readpdb","text":"readpdb(path) -> chains::Vector{ProteinChain}\n\n\n\n\n\n","category":"method"},{"location":"#ProteinChains.readrecord-Tuple{AbstractString, Type{<:Union{BioStructures.MMCIFFormat, BioStructures.PDBFormat}}}","page":"Home","title":"ProteinChains.readrecord","text":"readrecord(path, format) -> chains::ProteinStructure\n\nLoads a protein structure from a PDB file.\n\nExported formats: PDBFormat, MMCIFFormat\n\nExamples\n\nreadrecord(\"example.pdb\"); # detects PDB format from extension\n\nreadrecord(\"example.cif\"); # detects mmCIF format from extension\n\nreadrecord(\"example.abc\", PDBFormat); # force PDB format\n\nreadrecord(\"example.xyz\", MMCIFFormat); # force mmCIF format\n\n\n\n\n\n","category":"method"},{"location":"#ProteinChains.writechains-Tuple{AbstractString, AbstractVector{ProteinChain}, Type{BioStructures.PDBFormat}}","page":"Home","title":"ProteinChains.writechains","text":"writechains(path, chains::AbstractVector{ProteinChain}, format)\nwritechains(path, chain::ProteinChain, format)\n\nWrite a protein structure (represented as a Vector{ProteinChain}s) to file with the specified format.\n\nExported formats: PDBFormat, MMCIFFormat\n\nExamples\n\njulia> writechains(\"example.pdb\", chains) # detects PDB format from extension\n\njulia> writechains(\"example.cif\", chains) # detects mmCIF format from extension\n\njulia> writechains(\"example.abc\", chains, PDBFormat) # force PDB format\n\njulia> writechains(\"example.xyz\", chains, MMCIFFormat) # force mmCIF format\n\n\n\n\n\n","category":"method"},{"location":"#ProteinChains.writepdb-Tuple{AbstractString, AbstractVector{ProteinChain}}","page":"Home","title":"ProteinChains.writepdb","text":"writepdb(path, chains::AbstractVector{ProteinChain})\nwritepdb(path, chain::ProteinChain)\n\n\n\n\n\n","category":"method"}]
}
