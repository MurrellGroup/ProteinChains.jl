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
    print(io, summary(chain), ":")
    printstyled(io, "\n  fields:", color=:yellow)
    context = IOContext(io, :compact => true, :limit => true)
    for fieldname in fieldnames(ProteinChain)[1:end-1]
        printfield(context, chain, fieldname)
    end
    printstyled(io, "\n  properties:", color=:yellow)
    isempty(chain.properties) && print(io, " (none)")
    for property in keys(chain.properties)
        printfield(context, chain, property)
    end
end
