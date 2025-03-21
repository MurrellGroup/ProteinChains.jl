using Base.Filesystem

const MAX_CACHE_ENTRIES = 4
const CACHE_DIR = Ref{Union{String,Nothing}}(nothing)

function initialize_cache_dir()
    if isnothing(CACHE_DIR[])
        CACHE_DIR[] = mktempdir(prefix="pdb_cache_")
        atexit(() -> rm(CACHE_DIR[], recursive=true))
    end
end

function manage_cache()
    entries = readdir(CACHE_DIR[])
    while length(entries) > MAX_CACHE_ENTRIES
        oldest_entry = entries[argmin([mtime(joinpath(CACHE_DIR[], entry)) for entry in entries])]
        rm(joinpath(CACHE_DIR[], oldest_entry))
        entries = readdir(CACHE_DIR[])
    end
end

"""
    pdbentry(pdbid::AbstractString; format=MMCIFFormat, kws...)

Keyword arguments get propagated to [`BioStructures.downloadpdb`](https://biojulia.dev/BioStructures.jl/stable/api/#BioStructures.downloadpdb-Tuple{AbstractString})

Downloads are cached in a temporary directory.

## Examples

```jldoctest
julia> pdbentry("1EYE")
[ Info: Downloading file from PDB: 1EYE
1-element ProteinStructure "1EYE.cif":
 256-residue ProteinChain{Float64} (A)

julia> pdb"1EYE" # string macro for convenience
[ Info: File exists: 1EYE
1-element ProteinStructure "1EYE.cif":
 256-residue ProteinChain{Float64} (A)

julia> pdb"1EYE"A # string suffix to get a specific chain
[ Info: File exists: 1EYE
256-residue ProteinChain{Float64} (A)

julia> pdb"1EYE"1 # integer suffix to specify "ba_number" keyword
[ Info: Downloading file from PDB: 1EYE
2-element ProteinStructure "1EYE_ba1.cif":
 256-residue ProteinChain{Float64} (A)
 256-residue ProteinChain{Float64} (A-2)
```
"""
function pdbentry(pdbid::AbstractString; dir=nothing, format=MMCIFFormat, kws...)
    isnothing(dir) && return cachedpdbentry(pdbid; format, kws...)
    path = BioStructures.downloadpdb(pdbid; dir, format, kws...)
    return read(path, ProteinStructure, format)
end

function cachedpdbentry(pdbid::AbstractString; format=MMCIFFormat, kws...)
    initialize_cache_dir()
    structure = pdbentry(pdbid; dir=CACHE_DIR[], format, kws...)
    manage_cache()
    return structure
end

macro pdb_str(pdbid::AbstractString)
    :(pdbentry($(esc(pdbid))))
end

macro pdb_str(pdbid::AbstractString, chain::AbstractString)
    :(pdbentry($(esc(pdbid)))[$chain])
end

macro pdb_str(pdbid::AbstractString, ba_number::Integer)
    :(pdbentry($(esc(pdbid)), ba_number=$ba_number))
end

function getmmcifdict(pdbid::AbstractString; dir=nothing, kws...)
    isnothing(dir) && return cachedmmcifdict(pdbid; kws...)
    path = BioStructures.downloadpdb(pdbid; dir, format=MMCIFFormat, kws...)
    return MMCIFDict(path)
end

function cachedmmcifdict(pdbid::AbstractString; kws...)
    initialize_cache_dir()
    mmcifdict = getmmcifdict(pdbid; dir=CACHE_DIR[], kws...)
    manage_cache()
    return mmcifdict
end

macro mmcifdict_str(pdbid::AbstractString, ba_number::Integer=0)
    :(getmmcifdict($(esc(pdbid)), ba_number=$ba_number))
end
