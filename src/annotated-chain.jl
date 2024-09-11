@dynamic mutable struct AnnotatedProteinChain{T} <: AbstractProteinChain{T}
    id::String
    sequence::String
    backbone::Array{T,3}
    numbering::Vector{Int64}
    atoms::Vector{Vector{Atom{T}}}
    indexable_properties::Vector{Symbol}

    function AnnotatedProteinChain(
        id::AbstractString, sequence::AbstractString, backbone::AbstractArray{T,3},
        numbering::AbstractVector{<:Integer}, atoms::Vector{Vector{Atom{T}}},
        indexable_properties=Symbol[]; kwargs...
    ) where T
        @assert length(sequence) == size(backbone, 3) == length(atoms) == length(numbering)
        AnnotatedProteinChain{T}(id, sequence, backbone, atoms, numbering, indexable_properties; kwargs...)
    end
end

ProteinChain(chain::AbstractProteinChain) = ProteinChain(chain.id, chain.sequence, chain.backbone, chain.numbering, chain.atoms)
AnnotatedProteinChain(chain::AbstractProteinChain) = AnnotatedProteinChain(chain.id, chain.sequence, chain.backbone, chain.numbering, chain.atoms)

function Base.getindex(chain::AnnotatedProteinChain, i::AbstractVector{<:Integer})
    properties = Dict{Symbol,Any}()
    for name in getproperties(chain; fields=false)
        value = getproperty(chain, name)
        properties[name] = name in chain.indexable_properties ? collect(selectdim(value, ndims(value), i)) : value
    end
    return annotate(ProteinChain(chain)[i]; properties...)
end

function annotate_indexable!(chain::AnnotatedProteinChain; annotations...)
    for (name, value) in annotations
        @assert value isa AbstractArray
        @assert size(value)[end] == length(chain)
        push!(chain.indexable_properties, name)
        setproperty!(chain, name, value)
    end
    return chain
end

annotate_indexable(chain::AnnotatedProteinChain; annotations...) = annotate_indexable!(deepcopy(chain); annotations...)
annotate_indexable(chain::ProteinChain; annotations...) = annotate_indexable!(AnnotatedProteinChain(chain); annotations...)

function annotate!(chain::AnnotatedProteinChain; annotations...)
    for (name, value) in annotations
        setproperty!(chain, name, value)
    end
    return chain
end

annotate(chain::AnnotatedProteinChain; annotations...) = annotate(deepcopy(chain); annotations...)
annotate(chain::ProteinChain; annotations...) = annotate!(AnnotatedProteinChain(chain); annotations...)
