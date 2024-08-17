"""
    ProteinChain

## Examples
```jldoctest
julia> structure = pdb"1EYE"
[ Info: Downloading file from PDB: 1EYE
1-chain ProteinStructure "1EYE.cif" with 0 properties:
  256-residue ProteinChain "A" with 0 properties

julia> structure[1]
256-residue ProteinChain "A" with 0 properties:
  fields:
    id::String = "A"
    aminoacids::String = "PVQVMGVLNVTDDSFSDGGCYLDLDDAVKHGLAMAAAGAGIVDVGGETSRVIPVVKELAAQGITVSIDTMRADVARAALQNGAQMVNDVSGGRADPAM…
    backbone::Array{Float64, 3} = [45.592 44.171 43.719; -10.864 -10.936 -9.688; 30.192 30.504 31.278;;; 42.568 42.02 40.707; -9.163 …
    numbers::Vector{Int64} = [5, 6, 7, 8, 9, 10, 11, 12, 13, 14  …  265, 266, 267, 268, 269, 270, 271, 272, 273, 274]
  properties: (none)
```
"""
mutable struct ProteinChain
    id::String
    aminoacids::String
    backbone::Array{Float64,3}
    numbers::Vector{Int}
    properties::Dict{Symbol,Any}

    function ProteinChain(id::String, aminoacids::String, backbone::Array{Float64,3}, numbers::Vector{Int}, properties::Dict{Symbol,Any})
        chain = new(id, aminoacids, backbone, numbers, properties)
        countresidues(chain)
        @assert all(!in(fieldnames(ProteinChain)), keys(properties))
        return chain
    end
end

Base.getproperty(chain::ProteinChain, property::Symbol) = _getproperty(chain, property)
Base.setproperty!(chain::ProteinChain, property::Symbol, value) = _setproperty!(chain, property, value)

function countresidues(chain::ProteinChain)
    @assert length(chain.aminoacids) == size(chain.backbone, 3) == length(chain.numbers)
    return length(chain.aminoacids)
end

function ProteinChain(id::String, aminoacids::String, backbone::Array{Float64,3}, numbers::Vector{Int}; kwargs...)
    chain = ProteinChain(id, aminoacids, backbone, numbers, Dict{Symbol,Any}())
    !isempty(kwargs) && push!(chain.properties, kwargs...)
    return chain
end

Base.summary(chain::ProteinChain) = "$(countresidues(chain))-residue ProteinChain \"$(chain.id)\" with $(length(chain.properties)) properties"

function Base.show(io::IO, ::MIME"text/plain", chain::ProteinChain)
    context = IOContext(io, :compact => true, :limit => true)
    print(context, summary(chain), ":")
    printstyled(context, "\n  fields:", color=:yellow)
    for fieldname in fieldnames(ProteinChain)[1:end-1]
        printfield(context, chain, fieldname)
    end
    printstyled(context, "\n  properties:", color=:yellow)
    isempty(chain.properties) && print(io, " (none)")
    for property in keys(chain.properties)
        printfield(context, chain, property)
    end
end
