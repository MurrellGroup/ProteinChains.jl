expand_mask(L, freeze::Bool) = trues(L) .= !freeze
expand_mask(L, freeze::AbstractVector{Bool}) = .!freeze
function expand_mask(L, freeze::Union{UnitRange, AbstractVector{Integer}})
    mask = trues(L)
    mask[freeze] .= false
    return mask
end
get_aa_mask(ch::ProteinChain) = hasproperty(ch, :aafreeze) ? expand_mask(length(ch.numbering), ch.aafreeze) : trues(length(ch.numbering))
get_pos_mask(ch::ProteinChain) = hasproperty(ch, :posfreeze) ? expand_mask(length(ch.numbering), ch.posfreeze) : trues(length(ch.numbering))
get_flatmasks(prot::ProteinStructure) = vcat([get_pos_mask(ch) for ch in prot.chains]...), vcat([get_aa_mask(ch) for ch in prot.chains]...)
function masked_out_structure(prot::ProteinStructure{T}, flat_mask) where T
    latest = 0
    chains = ProteinChain{T}[]
    for ch in prot.chains
        mask_region = flat_mask[latest+1:latest+length(ch.numbering)]
        latest += length(ch.numbering)
        if sum(.!mask_region) > 0
            push!(chains, ch[findall(.!mask_region)])
        end
    end
    return ProteinStructure(prot.name, Atom{T}[], chains)
end
