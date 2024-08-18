# TODO: store non-backbone atoms

# TODO: property tree dependency structure

assign_property!(x, name::Symbol, args...) = assign_property!(x, Val(name), args...)

function assign_property!(chain::ProteinChain, ::Val{:secondary_structure}, dict::Union{Dict{ProteinChain,Any},Nothing}=nothing)
    number_vector = isnothing(dict) ? ASS.assign_secondary_structure(chain) : dict[chain]
    chain.secondary_structure = number_vector_to_code_string(number_vector)
end

function assign_property!(chain::ProteinChain, ::Val{:ideal_residue}, ideal_residue::AbstractMatrix{<:Real}=DEFAULT_IDEAL_RESIDUE)
    chain.ideal_residue = collect(ideal_residue)
end

function assign_property!(chain::ProteinChain, ::Val{:residue_rotations})
    frames = Backboner.Frames(Backboner.Backbone(chain.backbone), chain.ideal_residue)
    chain.residue_rotations = frames.rotations
end

function assign_property!(chain::ProteinChain, ::Val{:residue_rotations_quat})
    frames = Backboner.Frames(Backboner.Backbone(chain.backbone), chain.ideal_residue)
    chain.residue_rotations_quat = rotation_matrices_to_quaternions(frames.rotations)
end

function assign_property!(chain::ProteinChain, ::Val{:residue_translations})
    chain.residue_translations = Backboner.centroid(chain.backbone; dims=2)
end

function assign_property!(chain::ProteinChain, ::Val{:bond_lengths})
    chain.bond_lengths = Backboner.get_bond_lengths(Backboner.Backbone(chain.backbone))
end

function assign_property!(chain::ProteinChain, ::Val{:bond_angles})
    chain.bond_angles = Backboner.get_bond_angles(Backboner.Backbone(chain.backbone))
end

function assign_property!(chain::ProteinChain, ::Val{:torsion_angles})
    chain.torsion_angles = Backboner.get_torsional_angles(Backboner.Backbone(chain.backbone))
end

function assign_property!(chain::ProteinChain, ::Val{:is_knotted})
    chain.is_knotted = Backboner.is_knotted(Backboner.Backbone(chain.backbone)[2:3:end])
end


struct ResidueIndexingUndefined <: Exception
    name::Symbol
end

ResidueIndexingUndefined(val::Val) = typeof(val).type_parameters[1]

Base.showerror(io::IO, err::ResidueIndexingUndefined) = print(io, "$(typeof(err)): property `:$(err.name)` does not have a defined reorder function. Property needs reassignment.")

residue_indexing_function(name::Symbol) = residue_indexing_function(Val(name))

residue_indexing_function(val_name::Val) = throw(ResidueIndexingUndefined(val_name))

# move non-properties to getindex(::ProteinChain)
residue_indexing_function(::Val{:id}) = (x, i) -> x
residue_indexing_function(::Val{:sequence}) = (x, i) -> x[i]
residue_indexing_function(::Val{:backbone}) = (x, i) -> x[:, :, i]
residue_indexing_function(::Val{:atoms}) = (x, i) -> x[i]

residue_indexing_function(::Val{:secondary_structure}) = (x, i) -> x[i]
residue_indexing_function(::Val{:residue_rotations}) = (x, i) -> x[:, :, i]
residue_indexing_function(::Val{:residue_rotations_quat}) = (x, i) -> x[:, i]
residue_indexing_function(::Val{:residue_translations}) = (x, i) -> x[:, i]

flatten_property(chains::AbstractVector{<:ProteinChain}, ::Val{:sequence}) = join(chain -> chain.sequence, chains)


function extend_property! end