"""
    BackboneGeometry(;
        N_Ca_length = 1.46,
        Ca_C_length = 1.52,
        C_N_length = 1.33,

        N_Ca_C_angle = 1.94,
        Ca_C_N_angle = 2.03,
        C_N_Ca_angle = 2.13,
    )

Define the idealized bond lengths and bond angles of a protein backbone.
"""
@kwdef struct BackboneGeometry
    N_Ca_length = 1.46
    Ca_C_length = 1.52
    C_N_length = 1.33

    N_Ca_C_angle = 1.94
    Ca_C_N_angle = 2.04
    C_N_Ca_angle = 2.13
end

const DEFAULT_BACKBONE_GEOMETRY = BackboneGeometry()

"""
    IdealResidue{T<:AbstractFloat} <: AbstractMatrix{T}

    IdealResidue{T}(backbone_geometry=DEFAULT_BACKBONE_GEOMETRY; template=nothing) where T

A 3x3 matrix representing the idealized geometry of a protein residue, with columns representing
the N, Ca, and C atom positions of a residue positioned at the origin.
"""
struct IdealResidue{T<:AbstractFloat} <: AbstractMatrix{T}
    backbone_geometry::BackboneGeometry
    N_Ca_C_coords::Matrix{T}
end

function IdealResidue{T}(backbone_geometry::BackboneGeometry=DEFAULT_BACKBONE_GEOMETRY; template=nothing) where T
    N_Ca_C_coords = Matrix{T}(undef, 3, 3)
    N_pos, Ca_pos, C_pos = eachcol(N_Ca_C_coords)
    Θ = backbone_geometry.N_Ca_C_angle - π/2
    N_pos .= [0, 0, 0]
    Ca_pos .= [backbone_geometry.N_Ca_length, 0, 0]
    C_pos .= Ca_pos + backbone_geometry.Ca_C_length * [sin(Θ), cos(Θ), 0]
    N_Ca_C_coords .-= Backboner.centroid(N_Ca_C_coords; dims=2)
    if !isnothing(template)
        wanted_orientation, current_offset, wanted_offset = Backboner.kabsch_algorithm(N_Ca_C_coords, template)
        N_Ca_C_coords .= wanted_orientation * (N_Ca_C_coords .- current_offset) .+ wanted_offset
    end
    IdealResidue{T}(backbone_geometry, N_Ca_C_coords)
end

Base.size(::IdealResidue) = (3,3)
Base.getindex(ideal_residue::IdealResidue, args...) = ideal_residue.N_Ca_C_coords[args...]

"""
    STANDARD_RESIDUE_TEMPLATE

This is a template of a "standard residue", with a very specific and
distinct shape, size, and orientation. which needs to be consistent if we want to
represent protein structures as sets of residue rotations and translations.

Thus, we can use this residue as a template for aligning other residues with very precise
geometry to it.

```jldocotest
julia> IdealResidue{Float64}(BackboneGeometry(N_Ca_C_angle = 1.93); template=ProteinChains.STANDARD_RESIDUE_TEMPLATE)
3×3 IdealResidue{Float64}:
 -1.06447   -0.199174   1.26364
  0.646303  -0.529648  -0.116655
  0.0        0.0        0.0
```
"""
const STANDARD_RESIDUE_TEMPLATE = [
    -1.066  -0.200   1.266;
     0.645  -0.527  -0.118;
     0.000   0.000   0.000;
] #  N       Ca      C

const STANDARD_RESIDUE = IdealResidue{Float64}(DEFAULT_BACKBONE_GEOMETRY; template=STANDARD_RESIDUE_TEMPLATE)

"""
    append_residue(Backbone::Backbone, torsion_angles::Vector{<:Real}; ideal::BackboneGeometry=DEFAULT_BACKBONE_GEOMETRY)

Create a new backbone by appending 3 new torsion angles (ψ, ω, ϕ) at the end, using bond lengths and bond angles specified in `BackboneGeometry`.
"""
function append_residue(backbone::Backbone, torsion_angles::Vector{<:Real}; ideal::BackboneGeometry=DEFAULT_BACKBONE_GEOMETRY)
    length(torsion_angles) == 3 || throw(ArgumentError("length of `torsion_angles` must be 3, as only 3 backbone atoms are introduced."))
    bond_lengths = [ideal.C_N_length, ideal.N_Ca_length, ideal.Ca_C_length]
    bond_angles = [ideal.Ca_C_N_angle, ideal.C_N_Ca_angle, ideal.N_Ca_C_angle]
    append_bonds(backbone, Float64.(bond_lengths), Float64.(bond_angles), Float64.(torsion_angles))
end

"""
    append_residue(Backbone::Backbone, torsion_angles::Vector{<:Real}; ideal::BackboneGeometry=DEFAULT_BACKBONE_GEOMETRY)

Create a new backbone by prepending 3 new torsion angles (ψ, ω, ϕ) at the beginning, using bond lengths and bond angles specified in the `BackboneGeometry`.

!!! note
    The torsion angle order is the same as it would be when appending. The order is *not* reversed.
"""
function prepend_residue(backbone::Backbone, torsion_angles::Vector{<:Real}; ideal::BackboneGeometry=DEFAULT_BACKBONE_GEOMETRY)
    length(torsion_angles) == 3 || throw(ArgumentError("length of `torsion_angles` must be 3, as only 3 backbone atoms are introduced."))
    bond_lengths = [ideal.N_Ca_length, ideal.Ca_C_length, ideal.C_N_length]
    bond_angles = [ideal.N_Ca_C_angle, ideal.Ca_C_N_angle, ideal.C_N_Ca_angle]
    return prepend_bonds(backbone, Float64.(bond_lengths), Float64.(bond_angles), Float64.(torsion_angles))
end
