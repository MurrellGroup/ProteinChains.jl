KDTree = AssigningSecondaryStructure.KDTree
inrange = AssigningSecondaryStructure.inrange

"""
    all_atom_coords(chain)

Return a `3 × N` matrix of all-atom Cartesian coordinates for a `ProteinChain`.
This is written to be type-stable and allocation-friendly.
"""
function all_atom_coords(chain)
    n_atoms = sum(length, chain.atoms)
    n_atoms == 0 && return zeros(Float64, 3, 0)

    first_atom = chain.atoms[1][1]
    T = typeof(first_atom.x)
    coords = Array{T}(undef, 3, n_atoms)

    k = 1
    for residue_atoms in chain.atoms
        for atom in residue_atoms
            c = ProteinChains.atom_coords(atom)
            coords[1, k] = c[1]
            coords[2, k] = c[2]
            coords[3, k] = c[3]
            k += 1
        end
    end

    coords
end

"""
    structure_atom_coords(structure)

Flatten all atoms from all chains in `structure` into a single coordinate matrix
and parallel vectors giving the chain index and residue index for each atom.
"""
function structure_atom_coords(structure)
    n_atoms = 0
    for chain in structure
        n_atoms += sum(length, chain.atoms)
    end
    n_atoms == 0 && return zeros(Float64, 3, 0), Int[], Int[]

    # Assume all coordinates share a common numeric type
    first_atom = structure[1].atoms[1][1]
    T = typeof(first_atom.x)
    coords = Array{T}(undef, 3, n_atoms)
    atom_chain_idx = Vector{Int}(undef, n_atoms)
    atom_res_idx = Vector{Int}(undef, n_atoms)

    k = 1
    for (ci, chain) in enumerate(structure)
        for (ri, residue_atoms) in enumerate(chain.atoms)
            for atom in residue_atoms
                c = ProteinChains.atom_coords(atom)
                coords[1, k] = c[1]
                coords[2, k] = c[2]
                coords[3, k] = c[3]
                atom_chain_idx[k] = ci
                atom_res_idx[k] = ri
                k += 1
            end
        end
    end

    coords, atom_chain_idx, atom_res_idx
end

# === 1) Shape / compactness features ===

"""
    shape_features(chain)

For a single `ProteinChain`, return a named tuple with:
- `radius_of_gyration`: standard R_g based on all atoms;
- `mean_pairwise_distance`: mean distance between all atom pairs;
- `max_pairwise_distance`: maximum distance between any two atoms.
"""
function shape_features(chain::ProteinChain)
    coords = all_atom_coords(chain)
    n_atoms = size(coords, 2)
    if n_atoms <= 1
        return (radius_of_gyration = 0.0,
                mean_pairwise_distance = 0.0,
                max_pairwise_distance = 0.0)
    end

    # Center of mass and radius of gyration without allocations
    cx = zero(eltype(coords))
    cy = zero(eltype(coords))
    cz = zero(eltype(coords))
    for j in 1:n_atoms
        cx += coords[1, j]
        cy += coords[2, j]
        cz += coords[3, j]
    end
    inv_n = inv(eltype(coords)(n_atoms))
    cx *= inv_n
    cy *= inv_n
    cz *= inv_n

    rg2 = zero(eltype(coords))
    for j in 1:n_atoms
        dx = coords[1, j] - cx
        dy = coords[2, j] - cy
        dz = coords[3, j] - cz
        rg2 += dx*dx + dy*dy + dz*dz
    end
    rg = sqrt(rg2 * inv_n)

    # Pairwise distances (allocation-free)
    sum_d = 0.0
    count = 0
    max_d = 0.0
    for i in 1:(n_atoms - 1)
        xi = coords[1, i]
        yi = coords[2, i]
        zi = coords[3, i]
        for j in (i + 1):n_atoms
            dx = xi - coords[1, j]
            dy = yi - coords[2, j]
            dz = zi - coords[3, j]
            d = sqrt(dx*dx + dy*dy + dz*dz)
            sum_d += d
            count += 1
            d > max_d && (max_d = d)
        end
    end
    mean_d = count == 0 ? 0.0 : sum_d / count

    (radius_of_gyration = rg,
     mean_pairwise_distance = mean_d,
     max_pairwise_distance = max_d)
end

shape_features(structure) = [shape_features(chain) for chain in structure]

# === 2) Cysteine nature ===

"""
    cysteine_features(structure; distance_cutoff = 2.5)

For each chain, count:
- `num_paired::Int`: number of cysteines that participate in a disulfide bond
  (including those paired with cysteines in other chains),
- `num_cys::Int`: total number of cysteines in the chain.

Pairing is defined by SG–SG distances < `distance_cutoff` Å,
including pairs between different chains.
Returns a vector of named tuples, one per chain.
"""
function cysteine_features(structure; distance_cutoff = 2.5)
    sg_coords = Vector{Any}()
    sg_chain_idx = Int[]
    sg_res_idx = Int[]

    # Collect all cysteine SG atoms across the structure
    for (ci, chain) in enumerate(structure)
        for (ri, residue_atoms) in enumerate(chain.atoms)
            if chain.sequence[ri] != 'C'
                continue
            end
            sg_pos = findfirst(a -> ProteinChains.atom_name(a) == " SG ", residue_atoms)
            sg_pos === nothing && continue
            atom = residue_atoms[sg_pos]
            push!(sg_coords, ProteinChains.atom_coords(atom))
            push!(sg_chain_idx, ci)
            push!(sg_res_idx, ri)
        end
    end

    n_sg = length(sg_coords)
    paired_flags = falses(n_sg)
    if n_sg > 0
        coords_mat = reduce(hcat, sg_coords)
        tree = KDTree(coords_mat)
        for i in 1:n_sg
            idxs = inrange(tree, coords_mat[:, i], distance_cutoff)
            for j in idxs
                j == i && continue
                paired_flags[i] = true
                paired_flags[j] = true
            end
        end
    end

    # Map (chain, residue) -> paired?
    paired_map = Dict{Tuple{Int,Int},Bool}()
    for k in 1:n_sg
        paired_map[(sg_chain_idx[k], sg_res_idx[k])] = paired_flags[k]
    end

    # Build per-chain summary
    features = NamedTuple[]
    for (ci, chain) in enumerate(structure)
        cys_positions = findall(==( 'C'), collect(chain.sequence))
        num_cys = length(cys_positions)
        num_paired = num_cys == 0 ? 0 :
            count(pos -> get(paired_map, (ci, pos), false), cys_positions)
        push!(features, (num_paired = num_paired, num_cys = num_cys))
    end
    features
end

# === 3) Contact order ===

"""
    residue_pair_min_distances(chain; max_distance; use_numbering = false)

Compute the minimum all-atom distance for each residue pair in `chain` whose
minimum inter-atomic distance is ≤ `max_distance`, using a KDTree over all
atoms. Returns a `Dict{Tuple{Int,Int},Float64}` mapping an ordered residue
identifier pair `(i, j)` with `i < j` to the minimum distance in Å.

If `use_numbering == false` (default), the residue identifiers are 1-based
indices into the chain (i.e. positions in `chain.atoms`). If
`use_numbering == true`, the identifiers are taken from `chain.numbering`.
"""
function residue_pair_min_distances(chain::ProteinChain; max_distance, use_numbering = false)
    L = length(chain)
    numbering = chain.numbering
    # All-atom coordinates and atom-to-residue mapping
    coords = all_atom_coords(chain)
    if size(coords, 2) == 0
        return Dict{Tuple{Int,Int},Float64}()
    end

    atom_res_idx = Int[]
    for (ri, residue_atoms) in enumerate(chain.atoms)
        for _ in residue_atoms
            push!(atom_res_idx, ri)
        end
    end

    # Build KDTree once over all atoms
    tree = KDTree(coords)

    # Single inrange call, exactly as in the user's example.
    idxs_per_atom = inrange(tree, coords, max_distance)

    # Track minimum distance per residue-residue pair
    pair_dists = Dict{Tuple{Int,Int},Float64}()
    for i in 1:length(idxs_per_atom)
        res_i = atom_res_idx[i]
        id_i = use_numbering ? numbering[res_i] : res_i

        for j in idxs_per_atom[i]
            j <= i && continue
            res_j = atom_res_idx[j]
            res_i == res_j && continue

            id_j = use_numbering ? numbering[res_j] : res_j
            a = id_i < id_j ? id_i : id_j
            b = id_i < id_j ? id_j : id_i

            dx = coords[1, i] - coords[1, j]
            dy = coords[2, i] - coords[2, j]
            dz = coords[3, i] - coords[3, j]
            d = sqrt(dx*dx + dy*dy + dz*dz)

            key = (a, b)
            if haskey(pair_dists, key)
                d < pair_dists[key] && (pair_dists[key] = d)
            else
                pair_dists[key] = d
            end
        end
    end

    pair_dists
end

"""
    residue_contact_order(chain; distance_thresholds = (6.0, 8.0, 10.0))

Compute length-normalized contact order for a single chain, using all-atom contacts.

For each distance threshold `d`, contacts are residue pairs with any inter-atomic
distance < `d` Å, excluding pairs with sequence-number separation ≤ 1 based on
`chain.numbering`.

Contact order is defined as

  CO = (1 / (N * L_eff)) * sum_{(i,j) in C} |numbering[j] - numbering[i]|

where `L_eff` is `maximum(numbering) - minimum(numbering) + 1` and `N` is the
number of contacts at that threshold.
"""
function residue_contact_order(chain::ProteinChain; distance_thresholds = (6.0, 10.0), min_separation = 4)
    numbering = chain.numbering
    L = length(chain)
    if L <= 2
        names = Symbol.("co_$(replace(string(d), '.' => '_'))" for d in distance_thresholds)
        return NamedTuple{Tuple(names)}(fill(0.0, length(names)))
    end

    # Build a single residue-pair distance dictionary using the largest threshold
    max_d = maximum(distance_thresholds)
    pair_dists = residue_pair_min_distances(chain; max_distance = max_d, use_numbering = true)

    if isempty(pair_dists)
        names = Symbol.("co_$(replace(string(d), '.' => '_'))" for d in distance_thresholds)
        return NamedTuple{Tuple(names)}(fill(0.0, length(names)))
    end

    # Prepare accumulators per threshold
    n_thr = length(distance_thresholds)
    sums = zeros(Float64, n_thr)
    counts = zeros(Int, n_thr)

    # For each contacting residue pair, contribute to all thresholds it satisfies
    for ((i, j), d) in pair_dists
        sep = abs(j - i)
        sep <= min_separation && continue
        for (ti, thr) in enumerate(distance_thresholds)
            d <= thr || continue
            sums[ti] += sep
            counts[ti] += 1
        end
    end

    L_eff = maximum(numbering) - minimum(numbering) + 1
    cos = Float64[]
    for ti in 1:n_thr
        if counts[ti] == 0
            push!(cos, 0.0)
        else
            push!(cos, sums[ti] / (L_eff * counts[ti]))
        end
    end

    names = Symbol.("co_$(replace(string(d), '.' => '_'))" for d in distance_thresholds)
    NamedTuple{Tuple(names)}(Tuple(cos))
end

residue_contact_order(structure; distance_thresholds = (6.0, 10.0), min_separation = 4) =
    [residue_contact_order(chain; distance_thresholds = distance_thresholds, min_separation = min_separation) for chain in structure]

# === 4) Inter-chain contact counts ===

"""
    interchain_contact_counts(structure; distance_cutoff = 5.0)

For each chain, count how many residues make contact with at least one
residue in any other chain, based on all-atom contacts with a distance
cutoff `distance_cutoff` (Å).

Returns a vector of integers, one per chain.
"""
function interchain_contact_counts(structure; distance_cutoff = 5.0)
    n_chains = length(structure)
    n_chains <= 1 && return fill(0, n_chains)

    coords, atom_chain_idx, atom_res_idx = structure_atom_coords(structure)
    n_atoms = size(coords, 2)
    n_atoms == 0 && return fill(0, n_chains)

    tree = KDTree(coords)
    r = distance_cutoff

    # One boolean mask per chain, indexed by residue position
    has_contact_res = [falses(length(chain)) for chain in structure]

    for i in 1:n_atoms
        ci = atom_chain_idx[i]
        ri = atom_res_idx[i]
        has_contact_res[ci][ri] && continue

        point = @view coords[:, i]
        idxs = inrange(tree, point, r)
        for j in idxs
            j == i && continue
            cj = atom_chain_idx[j]
            cj == ci && continue
            has_contact_res[ci][ri] = true
            break
        end
    end

    [count(mask) for mask in has_contact_res]
end

"""
    interchain_residue_contact_proportion(structure; distance_cutoff = 5.0)

Residue-level inter-chain contact proportion: for each chain, returns the
fraction of residues that have at least one atom within `distance_cutoff` Å
of any atom in another chain.
"""
function interchain_residue_contact_proportion(structure; distance_cutoff = 5.0)
    counts = interchain_contact_counts(structure; distance_cutoff = distance_cutoff)
    [length(structure[i]) == 0 ? 0.0 : counts[i] / length(structure[i]) for i in 1:length(structure)]
end

"""
    atom_contact_order(chain; cutoff = 6.0)

Atom-level contact order (ALCO) matching the definition in Eric Alm's
`contactOrder.pl` script.

For a single chain:

- consider all non-hydrogen atoms,
- for each contacting atom pair within `cutoff` Å, with positive sequence
  separation in `chain.numbering`, accumulate that separation;
- the absolute CO is `order / counts`;
- the relative CO is `absolute_CO / (max(numbering) - min(numbering) + 1)`.

Results should match Plaxco et al. 1998, Table 1.
"""
function atom_contact_order(chain::ProteinChain; cutoff = 6.0)
    numbering = chain.numbering

    # All-atom coordinates and residue index per atom, skipping hydrogens
    coords_list = Vector{SVector{3,Float64}}()
    atom_res_idx = Int[]
    for (ri, residue_atoms) in enumerate(chain.atoms)
        for atom in residue_atoms
            ProteinChains.atom_symbol(atom) == "H" && continue
            push!(coords_list, ProteinChains.atom_coords(atom))
            push!(atom_res_idx, ri)
        end
    end

    if isempty(coords_list)
        return (absolute = 0.0, relative = 0.0)
    end

    coords = reduce(hcat, coords_list)
    tree = KDTree(coords)
    idxs_per_atom = inrange(tree, coords, cutoff)

    total_sep = 0.0
    counts = 0
    n_atoms = length(idxs_per_atom)

    for i in 1:n_atoms
        res_i = atom_res_idx[i]
        num_i = numbering[res_i]
        for j in idxs_per_atom[i]
            j == i && continue
            res_j = atom_res_idx[j]
            num_j = numbering[res_j]
            seq_dist = num_i - num_j
            seq_dist > 0 || continue
            counts += 1
            total_sep += seq_dist
        end
    end

    if counts == 0
        return (absolute = 0.0, relative = 0.0)
    end

    abs_co = total_sep / counts
    L_eff = maximum(numbering) - minimum(numbering) + 1
    rel_co = abs_co / L_eff
    (absolute = abs_co, relative = rel_co)
end

atom_contact_order(structure; cutoff = 6.0) =
    [atom_contact_order(chain; cutoff = cutoff) for chain in structure]

"""
    interchain_atom_contact_proportion(structure; distance_cutoff = 5.0)

Atom-level inter-chain contact proportion: for each chain, returns the
fraction of atoms that are within `distance_cutoff` Å of at least one atom
in another chain.

Results should match Plaxco et al. 1998, Table 1.
"""
function interchain_atom_contact_proportion(structure; distance_cutoff = 5.0)
    n_chains = length(structure)
    n_chains <= 1 && return fill(0.0, n_chains)

    coords, atom_chain_idx, _ = structure_atom_coords(structure)
    n_atoms = size(coords, 2)
    n_atoms == 0 && return fill(0.0, n_chains)

    tree = KDTree(coords)
    idxs_per_atom = inrange(tree, coords, distance_cutoff)

    has_contact = falses(n_atoms)
    for i in 1:n_atoms
        ci = atom_chain_idx[i]
        for j in idxs_per_atom[i]
            j == i && continue
            atom_chain_idx[j] == ci && continue
            has_contact[i] = true
            break
        end
    end

    totals = zeros(Int, n_chains)
    contacts = zeros(Int, n_chains)
    for i in 1:n_atoms
        ci = atom_chain_idx[i]
        totals[ci] += 1
        has_contact[i] && (contacts[ci] += 1)
    end

    [totals[ci] == 0 ? 0.0 : contacts[ci] / totals[ci] for ci in 1:n_chains]
end

