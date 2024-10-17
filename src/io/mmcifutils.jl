export getmmcif
export mapmmcif

function map_first_occurrence(u, v)
    d = Dict{eltype(u),eltype(v)}()
    for (x, y) in zip(u, v)
        haskey(d, x) || (d[x] = y)
    end
    d
end

map_last_occurrence(u, v) = Dict(zip(u, v))

compose_map(d1, d2, fallback="?") = Dict(k => get(d2, v, fallback) for (k,v) in d1)

getmmcif(mmcifdict::AbstractDict{String,Vector{String}}, key::AbstractString) = get(mmcifdict, key, String[])

"""
    mapmmcif(mmcifdict, field1 => field2, field3 => field4, ...)

```jldoctest
julia> import BioStructures

julia> filename BioStructures.downloadpdb("3HFM", format=BioStructures.MMCIFFormat);
[ Info: Downloading file from PDB: 3HFM

julia> mmcifdict = BioStructures.MMCIFDict(filename);

julia> mapmmcif(mmcifdict,
           "_atom_site.auth_asym_id"   => "_atom_site.label_entity_id",
           "_entity_src_gen.entity_id" => "_entity_src_gen.pdbx_gene_src_ncbi_taxonomy_id")
Dict{String, String} with 3 entries:
  "Y" => "9031"
  "L" => "10090"
  "H" => "10090"
```
"""
mapmmcif(mmcifdict, pairs::Pair{String,String}...) =
    mapreduce(((from,to),) -> map_first_occurrence(getmmcif(mmcifdict, from), getmmcif(mmcifdict, to)), compose_map, pairs)

get_auth_asym_to_entity(mmcifdict) = mapmmcif(mmcifdict, "_atom_site.auth_asym_id" => "_atom_site.label_entity_id")

function get_auth_asym_to_taxid(mmcifdict)
    mapmmcif(mmcifdict,
        "_atom_site.auth_asym_id"   => "_atom_site.label_entity_id",
        "_entity_src_gen.entity_id" => "_entity_src_gen.pdbx_gene_src_ncbi_taxonomy_id")
end
