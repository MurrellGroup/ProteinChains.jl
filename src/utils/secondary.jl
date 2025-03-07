AssigningSecondaryStructure.assign_secondary_structure(chain::ProteinChain) =
    AssigningSecondaryStructure.assign_secondary_structure(get_backbone(chain))

function AssigningSecondaryStructure.assign_secondary_structure!(chain::ProteinChain)
    ss = assign_secondary_structure(chain)
    chain.ss = ss
    return chain
end

function AssigningSecondaryStructure.assign_secondary_structure!(structure::ProteinStructure)
    for chain in structure
        assign_secondary_structure!(chain)
    end
end
