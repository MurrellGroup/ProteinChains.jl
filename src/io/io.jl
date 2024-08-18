using BioStructures: BioStructures, PDBFormat, MMCIFFormat

const ProteinFileFormat = Union{PDBFormat, MMCIFFormat}
const AMINOACIDS = Set("ACDEFGHIKLMNPQRSTVWY")
const BACKBONE_ATOM_NAMES = ["N", "CA", "C"]
const BACKBONE_ATOM_SYMBOLS = ["N", "C", "C"]

const threeletter_aa_names = Dict{Char, String}([Char(v) => k for (k, v) in BioStructures.threeletter_to_aa])
threeletter_to_oneletter(threeletter::AbstractString) = Char(get(BioStructures.threeletter_to_aa, threeletter, 'X'))
oneletter_resname(residue::BioStructures.AbstractResidue) = threeletter_to_oneletter(BioStructures.resname(residue))

const pdbextension_to_format = Dict(ext => format for (format, ext) in BioStructures.pdbextension)

get_format(path::AbstractString) = get(pdbextension_to_format, lowercase(last(splitext(path))[2:end]), PDBFormat)

include("read.jl")
include("write.jl")
include("download.jl")