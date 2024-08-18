using BioStructures: BioStructures, PDBFormat, MMCIFFormat

const ProteinFileFormat = Union{PDBFormat, MMCIFFormat}
const AMINOACIDS = Set("ACDEFGHIKLMNPQRSTVWY")
const BACKBONE_ATOM_NAMES = ["N", "CA", "C"]

threeletter_to_oneletter(threeletter::AbstractString) = Char(get(BioStructures.threeletter_to_aa, threeletter, 'X'))
oneletter_resname(residue::BioStructures.AbstractResidue) = threeletter_to_oneletter(BioStructures.resname(residue))

backbone_atom_selector(atom::BioStructures.AbstractAtom) = BioStructures.atomnameselector(atom, BACKBONE_ATOM_NAMES)
backbone_residue_selector(residue::BioStructures.AbstractResidue) =
    oneletter_resname(residue) in AMINOACIDS &&
    BioStructures.countatoms(residue, backbone_atom_selector) == 3 &&
    BioStructures.standardselector(residue) &&
    !BioStructures.disorderselector(residue)

function get_atom(residue::BioStructures.AbstractResidue, atom_name::AbstractString)
    selector = atom -> BioStructures.atomnameselector(atom, [atom_name])
    residue_atoms = BioStructures.collectatoms(residue, selector)
    return only(residue_atoms)
end

function get_backbone(residue::BioStructures.AbstractResidue)
    atom_name_to_atom_coords = atom_name -> BioStructures.coords(get_atom(residue, atom_name))
    residue_coords = stack(atom_name_to_atom_coords, BACKBONE_ATOM_NAMES)
    return Backbone(residue_coords)
end

function get_backbone(residues::Vector{BioStructures.AbstractResidue})
    chain_coords = mapreduce(residue -> get_backbone(residue).coords, hcat, residues; init=Matrix{Float64}(undef, 3, 0))
    return Backbone(chain_coords)
end

get_backbone(chain::BioStructures.Chain, selector) = get_backbone(BioStructures.collectresidues(chain, selector))

aminoacid_sequence(residues::Vector{BioStructures.AbstractResidue}) = join(oneletter_resname.(residues))
aminoacid_sequence(chain::BioStructures.Chain, selector) = aminoacid_sequence(BioStructures.collectresidues(chain, selector))

#=function get_residue_atoms(residues::Vector{BioStructures.AbstractResidue})
    atom_names = Vector{Vector{String}}[]
    atom_coords = Vector{Vector{Float64}}[]
    for residue in residues
        residue_atoms = BioStructures.collectatoms(residue, a -> !backbone_atom_selector(a) && BioStructures.standardselector(a) && !BioStructures.disorderselector(a))
        push!(atom_names, map(a -> a.name), residue_atoms)
        push!(atom_coords, map(a -> a.coords), residue_atoms)
    end
    return atom_names, atom_coords
end=#

function ProteinChain(residues::Vector{BioStructures.AbstractResidue})
    id = only(unique(BioStructures.chainid.(residues)))
    aminoacids = aminoacid_sequence(residues)
    backbone = reshape(Backbone(residues).coords, 3, 3, :)
    numbers = BioStructures.resnumber.(residues)
    #ins_codes = BioStructures.inscode.(residues)
    #modelnum = only(unique(BioStructures.modelnumber.(residues)))
    #residue_atoms = get_residue_atoms(residues)
    return ProteinChain(id, aminoacids, backbone, numbers)
end

function ProteinChain(chain::BioStructures.Chain, selector=backbone_residue_selector)
    residues = BioStructures.collectresidues(chain, selector)
    isempty(residues) && return ProteinChain(BioStructures.chainid(chain), "", zeros(3, 3, 0), Int[])
    return ProteinChain(residues)
end

function collectchains(struc::BioStructures.MolecularStructure, selector=backbone_residue_selector)
    chains = ProteinChain[]
    for model in struc, chain in model
        if !isempty(chain)
            pc = ProteinChain(chain, selector)
            countresidues(pc) > 0 && push!(chains, pc)
        end
    end
    return chains
end

function ProteinStructure(struc::BioStructures.MolecularStructure)
    name = struc.name
    chains = collectchains(struc)
    ProteinStructure(name, chains, Dict())
end

const pdbextension_to_format = Dict(ext => format for (format, ext) in BioStructures.pdbextension)

get_format(path::AbstractString) = get(pdbextension_to_format, lowercase(last(splitext(path))[2:end]), PDBFormat)

"""
    readchains(path, format) -> chains::ProteinStructure

Loads a protein structure from a PDB file.

Exported formats: `PDBFormat`, `MMCIFFormat`

## Examples

```julia
readchains("example.pdb") # detects PDB format from extension

readchains("example.cif") # detects mmCIF format from extension

readchains("example.abc", PDBFormat) # force PDB format

readchains("example.xyz", MMCIFFormat) # force mmCIF format
```
"""
readchains(path::AbstractString, format::Type{<:ProteinFileFormat}) = ProteinStructure(read(path, format))
readchains(path::AbstractString) = readchains(path, get_format(path))

"""
    readpdb(path) -> chains::Vector{ProteinChain}
"""
readpdb(path::AbstractString) = readchains(path, PDBFormat)

"""
    readcif(path) -> chains::Vector{ProteinChain}
"""
readcif(path::AbstractString) = readchains(path, MMCIFFormat)

"""
    writepdb(path, chains::AbstractVector{ProteinChain})
    writepdb(path, chain::ProteinChain)
"""
function writepdb(path::AbstractString, chains::AbstractVector{ProteinChain})
    atom_records = BioStructures.AtomRecord[]
    index = 0
    residue_index = 0
    for chain in chains
        residue_backbone_coords = reshape(chain.backbone.coords, 3, 3, :)
        assign_missing_oxygens!(chain)
        for i in 1:countresidues(chain)
            resname = get(threeletter_aa_names, chain.aminoacids[i], "XXX") # threletter_aa_names in residue.jl
            residue_index += 1
            residue_atoms = [
                [(; name, coords) for (name, coords) in zip(BACKBONE_ATOM_NAMES, eachcol(view(residue_backbone_coords, :, :, i)))];
                #filter(a -> !(a.name in BACKBONE_ATOM_NAMES), residue.atoms)
            ]
            for atom in residue_atoms
                index += 1
                number = chain.numbers[i]
                push!(atom_records, BioStructures.AtomRecord(false, index, atom.name, ' ', resname, chain.id,
                    number, ' ', atom.coords, 1.0, 0.0, strip(atom.name)[1:1], ""))
            end
        end
    end
    pdblines = BioStructures.pdbline.(atom_records)
    open(path, "w") do io
        for line in pdblines
            println(io, line)
        end
    end
    return nothing
end

writepdb(path::AbstractString, chain::ProteinChain) = writepdb(path, [chain])

# compat
writepdb(chains::Union{ProteinChain, AbstractVector{ProteinChain}}, path::AbstractString) = writepdb(path, chains)

"""
    writechains(path, chains::AbstractVector{ProteinChain}, format)
    writechains(path, chain::ProteinChain, format)

Write a protein structure (represented as a `Vector{ProteinChain}`s) to file with the specified format.

Exported formats: `PDBFormat`, `MMCIFFormat`

## Examples

```julia
writechains("example.pdb", chains) # detects PDB format from extension

writechains("example.cif", chains) # detects mmCIF format from extension

writechains("example.abc", chains, PDBFormat) # force PDB format

writechains("example.xyz", chains, MMCIFFormat) # force mmCIF format
```
"""
writechains(path::AbstractString, chains::AbstractVector{ProteinChain}, ::Type{PDBFormat}) = writepdb(path, chains)

function writechains(path::AbstractString, chains::AbstractVector{ProteinChain}, format::Type{<:ProteinFileFormat})
    struc = mktempdir() do temp_dir
        temp_path = joinpath(temp_dir, "temp.pdb")
        writechains(temp_path, chains, PDBFormat)
        read(temp_path, PDBFormat) # loads BioStructures.MolecularStructure
    end
    write_function = format == PDBFormat ? BioStructures.writepdb : BioStructures.writemmcif
    write_function(path, struc)
end

writechains(path::AbstractString, chains::AbstractVector{ProteinChain}) = writechains(path, chains, get_format(path))

writechains(path, chain::ProteinChain, args...) = writechains(path, [chain], args...)

pdbentry(pdbid::AbstractString; format=BioStructures.MMCIFFormat, kwargs...) = mktempdir() do dir
    path = BioStructures.downloadpdb(pdbid; dir, format, kwargs...)
    readchains(path, format)
end

macro pdb_str(pdbid)
    quote
        pdbentry($(esc(pdbid)))
    end
end