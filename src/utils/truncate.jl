using Statistics: mean
using LinearAlgebra: norm

function truncate(structure::ProteinStructure, residue_threshold::Integer)
    function _truncate(structure, selected=1:length(structure))
        new_structure = structure[selected]
        new_structure.truncated = length(selected) < length(structure)
        return new_structure
    end

    sum(length, structure; init=0) <= residue_threshold && return _truncate(structure)

    centers = [mean(get_backbone(chain), dims=(2,3)) for chain in structure]
    lengths = length.(structure)

    root = argmin(lengths)

    selected = Int[]
    total_residues = 0
    
    n_chains = length(structure)
    while total_residues < residue_threshold && length(selected) < n_chains
        distances = [norm(centers[i] - centers[root]) for i in 1:n_chains]
        distances[selected] .= Inf
        next_chain = argmin(distances)
        total_residues + lengths[next_chain] > residue_threshold && break
        push!(selected, next_chain)
        total_residues += lengths[next_chain]
    end

    return _truncate(structure, selected)
end
