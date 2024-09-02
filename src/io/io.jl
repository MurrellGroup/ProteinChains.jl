using BioStructures: BioStructures, PDBFormat, MMCIFFormat

const ProteinFileFormat = Union{PDBFormat, MMCIFFormat}
const AMINOACIDS = Set("ACDEFGHIKLMNPQRSTVWY")
const BACKBONE_ATOM_NAMES = ["N", "CA", "C"]
const BACKBONE_ATOM_SYMBOLS = ["N", "C", "C"]

const aa_to_threeletter = Dict{Char, String}([Char(v) => k for (k, v) in BioStructures.threeletter_to_aa])
threeletter_resname(aa::Char) = get(aa_to_threeletter, aa, "XXX")
oneletter_resname(threeletter::AbstractString) = Char(get(BioStructures.threeletter_to_aa, threeletter, 'X'))
oneletter_resname(residue::BioStructures.AbstractResidue) = oneletter_resname(BioStructures.resname(residue))

const pdbextension_to_format = Dict(ext => format for (format, ext) in BioStructures.pdbextension)

get_format(path::AbstractString) = get(pdbextension_to_format, lowercase(last(splitext(path))[2:end]), PDBFormat)

include("renumber.jl")
include("read.jl")
include("write.jl")

pdbentry(pdbid::AbstractString; format=BioStructures.MMCIFFormat, kwargs...) = mktempdir() do dir
    path = BioStructures.downloadpdb(pdbid; dir, format, kwargs...)
    readchains(path, format)
end

macro pdb_str(pdbid::AbstractString)
    :(pdbentry($(esc(pdbid))))
end
