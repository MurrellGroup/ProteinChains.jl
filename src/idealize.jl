Base.@kwdef struct BackboneGeometry
    N_Ca_length = 1.46 
    Ca_C_length = 1.52
    C_N_length = 1.33

    N_Ca_C_angle = 1.94
    Ca_C_N_angle = 2.03
    C_N_Ca_angle = 2.13
end

# TODO: make this align with default BackboneGeometry
const STANDARD_RESIDUE_ANGSTROM = [
    -1.066  -0.200   1.266;
     0.645  -0.527  -0.118;
     0.000   0.000   0.000;
] #  N       Ca      Cs

function append_residue(backbone::Backbone, dihedrals::Vector{<:Real}; ideal::BackboneGeometry=BackboneGeometry())
    length(dihedrals) == 3 || throw(ArgumentError("length of `dihedrals` must be 3, as only 3 backbone atoms are introduced."))
    lengths = [ideal.C_N_length, ideal.N_Ca_length, ideal.Ca_C_length]
    angles = [ideal.Ca_C_N_angle, ideal.C_N_Ca_angle, ideal.N_Ca_C_angle]
    append_bonds(backbone, Float64.(lengths), Float64.(angles), Float64.(dihedrals))
end

function prepend_residue(backbone::Backbone, dihedrals::Vector{<:Real}; ideal::BackboneGeometry=BackboneGeometry())
    length(dihedrals) == 3 || throw(ArgumentError("length of `dihedrals` must be 3, as only 3 backbone atoms are introduced."))
    lengths = [ideal.N_Ca_length, ideal.Ca_C_length, ideal.C_N_length]
    angles = [ideal.N_Ca_C_angle, ideal.Ca_C_N_angle, ideal.C_N_Ca_angle]
    return prepend_bonds(backbone, Float64.(lengths), Float64.(angles), Float64.(dihedrals))
end

# TODO: use cooldown (and Float64?) to increase precision and accuracy
# See ext/ZygoteExt.jl
function idealize! end
