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
Base.@kwdef struct BackboneGeometry
    N_Ca_length = 1.46
    Ca_C_length = 1.52
    C_N_length = 1.33

    N_Ca_C_angle = 1.94
    Ca_C_N_angle = 2.03
    C_N_Ca_angle = 2.13
end

const DEFAULT_BACKBONE_GEOMETRY = BackboneGeometry()

"""
    IdealResidue{T<:AbstractFloat} <: AbstractMatrix{T}

    IdealResidue{T}(backbone_geometry::BackboneGeometry=DEFAULT_BACKBONE_GEOMETRY)
"""
struct IdealResidue{T<:AbstractFloat} <: AbstractMatrix{T}
    backbone_geometry::BackboneGeometry
    N_Ca_C_coords::Matrix{T}

    function IdealResidue{T}(backbone_geometry::BackboneGeometry=DEFAULT_BACKBONE_GEOMETRY) where T
        N_Ca_C_coords = Matrix{T}(undef, 3, 3)
        N_pos, Ca_pos, C_pos = eachcol(N_Ca_C_coords)
        Θ = backbone_geometry.N_Ca_C_angle - π/2
        N_pos .= [0, 0, 0]
        Ca_pos .= [backbone_geometry.N_Ca_length, 0, 0]
        C_pos .= Ca_pos + backbone_geometry.Ca_C_length * [sin(Θ), cos(Θ), 0]
        centroid = (N_pos + Ca_pos + C_pos) / 3
        N_Ca_C_coords .-= centroid
        new{T}(backbone_geometry, N_Ca_C_coords)
    end
end

Base.size(::IdealResidue) = (3,3)
Base.getindex(ideal_residue::IdealResidue, args...) = ideal_residue.N_Ca_C_coords[args...]

const DEFAULT_IDEAL_RESIDUE = IdealResidue{Float64}()

const OLD_STANDARD_RESIDUE_ANGSTROM = [
    -1.066  -0.200   1.266;
     0.645  -0.527  -0.118;
     0.000   0.000   0.000;
] #  N       Ca      Cs

"""
    append_residue(Backbone::Backbone, torsion_angles::Vector{<:Real}; ideal::BackboneGeometry=DEFAULT_BACKBONE_GEOMETRY)

Create a new backbone by appending 3 torsion angles (ψ, ω, ϕ) at the end, using bond lengths and bond angles specified in BackboneGeometry.
"""
function append_residue(backbone::Backbone, torsion_angles::Vector{<:Real}; ideal::BackboneGeometry=DEFAULT_BACKBONE_GEOMETRY)
    length(torsion_angles) == 3 || throw(ArgumentError("length of `torsion_angles` must be 3, as only 3 backbone atoms are introduced."))
    bond_lengths = [ideal.C_N_length, ideal.N_Ca_length, ideal.Ca_C_length]
    bond_angles = [ideal.Ca_C_N_angle, ideal.C_N_Ca_angle, ideal.N_Ca_C_angle]
    append_bonds(backbone, Float64.(bond_lengths), Float64.(bond_angles), Float64.(torsion_angles))
end

"""
    append_residue(Backbone::Backbone, torsion_angles::Vector{<:Real}; ideal::BackboneGeometry=DEFAULT_BACKBONE_GEOMETRY)

Create a new backbone by prepending 3 torsion angles (ψ, ω, ϕ) at the end, using bond lengths and bond angles specified in BackboneGeometry.

!!! note
    The torsion angle order is the same as it would be when appending. The order is *not* reversed.
"""
function prepend_residue(backbone::Backbone, torsion_angles::Vector{<:Real}; ideal::BackboneGeometry=DEFAULT_BACKBONE_GEOMETRY)
    length(torsion_angles) == 3 || throw(ArgumentError("length of `torsion_angles` must be 3, as only 3 backbone atoms are introduced."))
    bond_lengths = [ideal.N_Ca_length, ideal.Ca_C_length, ideal.C_N_length]
    bond_angles = [ideal.N_Ca_C_angle, ideal.Ca_C_N_angle, ideal.C_N_Ca_angle]
    return prepend_bonds(backbone, Float64.(bond_lengths), Float64.(bond_angles), Float64.(torsion_angles))
end

# TODO: use cooldown (and Float64?) to increase precision and accuracy
# See ext/ZygoteExt.jl
function idealize! end
