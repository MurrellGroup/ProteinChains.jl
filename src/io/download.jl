using Base.Filesystem

const MAX_CACHE_ENTRIES = 4
const CACHE_DIR = Ref{Union{String,Nothing}}(nothing)

function initialize_cache_dir()
    if CACHE_DIR[] === nothing
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

@doc raw"""
    pdbentry(pdbid::AbstractString; format=MMCIFFormat, kwargs...)

Keyword arguments get propagated to [`BioStructures.downloadpdb`](https://biojulia.dev/BioStructures.jl/stable/api/#BioStructures.downloadpdb-Tuple{AbstractString})

Downloads are cached in a temporary directory up to $MAX_CACHE_ENTRIES times.

## Examples

```jldoctest
julia> pdbentry("1EYE")
[ Info: Downloading file from PDB: 1EYE
1-chain ProteinStructure{Float64} "1EYE.cif"
 256-residue ProteinChain{Float64, @NamedTuple{}} (A)

julia> pdb"1EYE" # convenience macro
[ Info: File exists: 1EYE
1-chain ProteinStructure{Float64} "1EYE.cif"
 256-residue ProteinChain{Float64, @NamedTuple{}} (A)

julia> pdb"1EYE"A # string suffix to get a specific chain
[ Info: File exists: 1EYE
256-residue ProteinChain{Float64, @NamedTuple{}} (A)

julia> pdb"1EYE"1 # integer suffix to specify "ba_number" keyword
[ Info: Downloading file from PDB: 1EYE
2-chain ProteinStructure{Float64} "1EYE_ba1.cif"
 256-residue ProteinChain{Float64, @NamedTuple{}} (A)
 256-residue ProteinChain{Float64, @NamedTuple{}} (A-2)
```
"""
function pdbentry(pdbid::AbstractString; format=MMCIFFormat, kwargs...)
    initialize_cache_dir()
    path = BioStructures.downloadpdb(pdbid; dir=CACHE_DIR[], format=format, kwargs...)
    manage_cache()
    return read(path, ProteinStructure, format)
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
