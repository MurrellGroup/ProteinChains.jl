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

function sequence_mismatch(master::ProteinChain, target::ProteinChain, master_inds, target_inds)
    matches = count(splat((m,t) -> master.sequence[m] == target.sequence[t]), zip(master_inds, target_inds); init=0)
    return 1 - matches / length(target_inds)
end

function detect_repeats(structure::ProteinStructure; mismatch_tolerance=0.05, length_threshold=0.8)
    N = length(structure.chains)
    L = sum(length, structure.chains)
    graph = SimpleGraph(L)
    assigned_chains = Set{Int}()
    offsets = [0; accumulate(+, length.(structure.chains))]
    for i in 1:N
        i in assigned_chains && continue
        master_chain = structure.chains[i]
        for j in i+1:N
            j in assigned_chains && continue
            target_chain = structure.chains[j]
            /(minmax(length(master_chain), length(target_chain))...) < length_threshold && continue
            master_inds, target_inds = residue_mapping(master_chain, target_chain)
            if sequence_mismatch(master_chain, target_chain, master_inds, target_inds) ≤ mismatch_tolerance
                for (t, m) in zip(target_inds, master_inds)
                    add_edge!(graph, offsets[i] + t, offsets[j] + m)
                end
                push!(assigned_chains, j)
            end
        end
    end
    return graph
end
