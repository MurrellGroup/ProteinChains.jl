import AssigningSecondaryStructure: assign_secondary_structure!, assign_secondary_structure

export assign_secondary_structure!

const SS_CODES = ['-', 'H', 'E']

number_vector_to_codes(v::AbstractVector{<:Integer}) = SS_CODES[v]
number_vector_to_code_string(v::AbstractVector{<:Integer}) = join(SS_CODES)[v]

function assign_secondary_structure!(chain::ProteinChain)
    number_vector = assign_secondary_structure(chain.backbone)
    chain.secondary_structure = number_vector_to_code_string(number_vector)
    return chain
end

function assign_secondary_structure!(structure::AbstractVector{<:ProteinChain})
    assign_secondary_structure!.(structure)
    return structure
end
