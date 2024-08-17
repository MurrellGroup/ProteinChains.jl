# TODO: store non-backbone atoms

function assign_bond_lengths!(chain::ProteinChain)
    chain.bond_lengths = Backboner.get_bond_lengths(Backboner.Backbone(chain.backbone))
end

function assign_bond_angles!(chain::ProteinChain)
    chain.bond_angles = Backboner.get_bond_angles(Backboner.Backbone(chain.backbone))
end

function assign_torsion_angles!(chain::ProteinChain)
    chain.torsion_angles = Backboner.get_torsional_angles(Backboner.Backbone(chain.backbone))
end

function assign_secondary_structure!(chain::ProteinChain, dict::Union{Dict{ProteinChain,Any},Nothing}=nothing)
    number_vector = isnothing(dict) ? ASS.assign_secondary_structure(chain) : dict[chain]
    chain.secondary_structure = number_vector_to_code_string(number_vector)
end

function assign_standard_residue!(chain::ProteinChain, standard_residue::Matrix{<:Real}=STANDARD_RESIDUE_ANGSTROM)
    chain.standard_residue = standard_residue
end

function assign_residue_rotations!(chain::ProteinChain)
    frames = Backboner.Frames(Backboner.Backbone(chain.backbone), chain.standard_residue)
    chain.residue_rotations = frames.rotations
end

function assign_residue_rotations_quat!(chain::ProteinChain)
    frames = Backboner.Frames(Backboner.Backbone(chain.backbone), chain.standard_residue)
    chain.residue_rotations_quat = rotation_matrices_to_quaternions(frames.rotations)
end

function assign_residue_translations!(chain::ProteinChain)
    chain.residue_translations = Backboner.centroid(chain.backbone; dims=2)
end

function assign_is_knotted!(chain::ProteinChain)
    chain.is_knotted = Backboner.is_knotted(Backboner.Backbone(chain.backbone)[2:3:end])
end
