using Base.Filesystem

const MAX_CACHE_ENTRIES = Ref{Int}(4)
const CACHE_DIR = Ref{Union{String, Nothing}}(nothing)

function initialize_cache_dir()
    if CACHE_DIR[] === nothing
        CACHE_DIR[] = mktempdir(prefix="pdb_cache_")
        atexit(() -> rm(CACHE_DIR[], recursive=true))
    end
end

function manage_cache()
    entries = readdir(CACHE_DIR[])
    while length(entries) > MAX_CACHE_ENTRIES[]
        oldest_entry = entries[argmin([mtime(joinpath(CACHE_DIR[], entry)) for entry in entries])]
        rm(joinpath(CACHE_DIR[], oldest_entry))
        entries = readdir(CACHE_DIR[])
    end
end

function pdbentry(pdbid::AbstractString; format=BioStructures.MMCIFFormat, kwargs...)
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
