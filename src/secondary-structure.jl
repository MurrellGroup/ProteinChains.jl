import AssigningSecondaryStructure as ASS

const SS_CODES = ['-', 'H', 'E']

number_vector_to_codes(v::AbstractVector{<:Integer}) = SS_CODES[v]
number_vector_to_code_string(v::AbstractVector{<:Integer}) = join(SS_CODES)[v]

function ASS.assign_secondary_structure(structure::AbstractVector{ProteinChain})
    chains = collect(structure)
    backbones = map(chain -> chain.backbone, chains)
    return Dict(zip(chains, ASS.assign_secondary_structure(backbones)))
end

ASS.assign_secondary_structure(chain::ProteinChain) = ASS.assign_secondary_structure([chain])[chain]