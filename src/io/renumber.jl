"""
    renumber!(structure::ProteinStructure, mmcif_dict::BioStructures.MMCIFDict)

Renumber the residues in a `ProteinStructure` object according to the numbering aligned to a reference sequence in the MMCIF file.
"""
function renumber!(structure::ProteinStructure, mmcif_dict::BioStructures.MMCIFDict)
    label_seq_ids = mmcif_dict["_atom_site.label_seq_id"]
    auth_seq_ids = mmcif_dict["_atom_site.auth_seq_id"]
    auth_asym_ids = mmcif_dict["_atom_site.auth_asym_id"]

    reduced = Int[]
    current_auth_seq_id = ""
    for (i, auth_seq_id) in enumerate(auth_seq_ids)
        current_auth_seq_id == auth_seq_id && continue
        current_auth_seq_id = auth_seq_id
        push!(reduced, i)
    end

    label_seq_ids = label_seq_ids[reduced]
    auth_seq_ids = auth_seq_ids[reduced]
    auth_asym_ids = auth_asym_ids[reduced]

    chainwise_auth_to_label_seq_ids = Dict{String,Dict{String,String}}()

    for asym_id in unique(auth_asym_ids)
        indices = auth_asym_ids .== asym_id
        chainwise_auth_to_label_seq_ids[asym_id] = Dict(zip(auth_seq_ids[indices], label_seq_ids[indices]))
    end

    for chain in structure
        renum_dict = chainwise_auth_to_label_seq_ids[chain.id]
        numbering_str = map(n -> get(renum_dict, n, "?"), string.(chain.numbering))
        "?" in numbering_str && break
        chain.numbering = parse.(Int, numbering_str)
    end

    return structure
end

#=function renumber!(structure::ProteinStructure, mmcif_dict::BioStructures.MMCIFDict)
    asym_ids = mmcif_dict["_pdbx_poly_seq_scheme.asym_id"]
    auth_seq_nums = map(s -> s == "?" ? -1 : parse(Int, s), mmcif_dict["_pdbx_poly_seq_scheme.auth_seq_num"])
    ndb_seq_nums = map(s -> s == "?" ? -1 : parse(Int, s), mmcif_dict["_pdbx_poly_seq_scheme.ndb_seq_num"])
    amino_acids = map(s -> get(BioStructures.threeletter_to_aa, s, BioStructures.AA_X), mmcif_dict["_pdbx_poly_seq_scheme.mon_id"])

    chainwise_auth_to_ndb = Dict{String,Dict{Int,Int}}()
    chainwise_sequences = Dict{String,String}()

    for asym_id in unique(asym_ids)
        indices = asym_ids .== asym_id
        chainwise_auth_to_ndb[asym_id] = Dict(zip(auth_seq_nums[indices], ndb_seq_nums[indices]))
        chainwise_sequences[asym_id] = join(amino_acids[indices])
    end

    for chain in structure
        chain.numbering = map(n -> chainwise_auth_to_ndb[chain.id][n], chain.numbering)
        chain.sequence_reference = chainwise_sequences[chain.id]
    end

    return structure
end=#