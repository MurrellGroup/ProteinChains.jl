import AssigningSecondaryStructure: assign_secondary_structure!, assign_secondary_structure

const SS_CODES = ['-', 'H', 'E']

number_vector_to_codes(v::AbstractVector{<:Integer}) = SS_CODES[v]

assign_secondary_structure(chain::AbstractProteinChain) = assign_secondary_structure(chain.backbone)

function assign_secondary_structure!(chain::AnnotatedProteinChain)
    number_vector = assign_secondary_structure(chain)
    annotate_indexable!(chain; secondary_structure=number_vector_to_codes(number_vector))
    return chain
end
