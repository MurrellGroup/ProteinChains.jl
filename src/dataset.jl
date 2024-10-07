"""
    ProteinDataset <: AbstractDict{String,ProteinStructure}
"""
struct ProteinDataset <: AbstractDict{String,ProteinStructure}
    dict::Dict{String,ProteinStructure}
end

ProteinDataset(structures::AbstractVector{<:ProteinStructure}) =
    ProteinDataset(Dict(structure.name => structure for structure in structures))

Base.length(d::ProteinDataset) = length(d.dict)
Base.iterate(d::ProteinDataset, args...) = iterate(d.dict, args...)
Base.:(==)(d1::ProteinDataset, d2::ProteinDataset) = d1.dict == d2.dict
