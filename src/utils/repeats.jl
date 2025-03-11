using Graphs

function writeproperty(lazy::Lazy{<:ProteinStructure}, ::Val{:repeats}, graph::SimpleGraph)
    repeats_group = HDF5.create_group(lazy.group, "repeats")
    repeats_group["nv"] = nv(graph)
    repeats_group["src"] = src.(edges(graph))
    repeats_group["dst"] = dst.(edges(graph))
end

function readproperty(lazy::Lazy, ::Val{:repeats})
    repeats_group = lazy.group["repeats"]
    nv = read(repeats_group, "nv")
    src = read(repeats_group, "src")
    dst = read(repeats_group, "dst")
    graph = SimpleGraph(nv)
    for (s, d) in zip(src, dst)
        add_edge!(graph, s, d)
    end
    return graph
end

function residue_mapping(master::ProteinChain, target::ProteinChain)
    master_map = Dict(master.numbering .=> 1:length(master.sequence))
    target_map = Dict(target.numbering .=> 1:length(target.sequence))
    common_resnums = intersect(keys(master_map), keys(target_map))
    return [master_map[i] for i in common_resnums], [target_map[i] for i in common_resnums]
end

mismatches(master::ProteinChain, target::ProteinChain, master_inds, target_inds) =
    count(splat((m,t) -> master.sequence[m] != target.sequence[t]), zip(master_inds, target_inds); init=0)

function detect_repeats(structure::ProteinStructure;
    mismatch_tolerance=1,
    length_threshold=0.8, # backwards compatibility
    overlap_proportion=length_threshold,
)
    N = length(structure)
    L = sum(length, structure)
    graph = SimpleGraph(L)
    assigned_chains = Set{Int}()
    offsets = [0; accumulate(+, length.(structure))]
    for i in 1:N
        i in assigned_chains && continue
        master_chain = structure[i]
        for j in i+1:N
            j in assigned_chains && continue
            target_chain = structure[j]
            splat(/)(minmax(length(master_chain), length(target_chain))) < overlap_proportion && continue
            master_inds, target_inds = residue_mapping(master_chain, target_chain)
            length(target_inds) / length(target_chain) < overlap_proportion && continue
            length(master_inds) / length(master_chain) < overlap_proportion && continue
            if mismatches(master_chain, target_chain, master_inds, target_inds) ≤ mismatch_tolerance
                for (t, m) in zip(target_inds, master_inds)
                    add_edge!(graph, offsets[i] + m, offsets[j] + t)
                end
                push!(assigned_chains, j)
            end
        end
    end
    return graph
end
