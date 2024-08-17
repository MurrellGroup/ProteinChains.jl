module ZygoteExt

using ProteinChains

import Zygote

function idealize!(chain::ProteinChain; ideal::BackboneGeometry=BackboneGeometry())
    backbone = Backboner.Backbone(chain.backbone)
    T = eltype(backbone.coords)
    ideal_lengths = T[ideal.N_Ca_length, ideal.Ca_C_length, ideal.C_N_length]
    ideal_angles = T[ideal.N_Ca_C_length, ideal.Ca_C_N_length, ideal.C_N_Ca_length]
    new_backbone = Backboner.idealize(backbone, ideal_lengths, ideal_angles)
    chain.backbone = new_backbone
    return chain
end

function idealize!(structure::ProteinStructure; kwargs...)
    for chain in structure
        idealize!(chain; kwargs...)
    end
    return structure
end

end
