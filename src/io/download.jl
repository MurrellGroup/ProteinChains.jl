pdbentry(pdbid::AbstractString; format=BioStructures.MMCIFFormat, kwargs...) = mktempdir() do dir
    path = BioStructures.downloadpdb(pdbid; dir, format, kwargs...)
    readchains(path, format)
end

macro pdb_str(pdbid)
    quote
        pdbentry($(esc(pdbid)))
    end
end
