const missingvals = Set([".", "?"])

function _renumber(structure::ProteinStructure, mmcif_dict::BioStructures.MMCIFDict)
    label_seq_ids = mmcif_dict["_atom_site.label_seq_id"]
    pdbx_PDB_ins_codes = mmcif_dict["_atom_site.pdbx_PDB_ins_code"]
    auth_seq_ids = mmcif_dict["_atom_site.auth_seq_id"]
    auth_asym_ids = mmcif_dict["_atom_site.auth_asym_id"]

    label_seq_id_dict = Dict{String,String}()

    for asym_id in unique(auth_asym_ids)
        indices = findall(auth_asym_ids .== asym_id)
        for i in indices
            ins_code = pdbx_PDB_ins_codes[i] == "?" ? " " : pdbx_PDB_ins_codes[i]
            key = asym_id*auth_seq_ids[i]*ins_code
            !haskey(label_seq_id_dict, key) && (label_seq_id_dict[key] = label_seq_ids[i])
        end
    end

    chainwise_renumbering = Vector{Int32}[]
    for chain in structure
        if length(chain) == 0
            push!(chainwise_renumbering, Int32[])
            continue
        end
        renumbering_str = map(
            (resnum, ins_code) -> parse(Int32, get(label_seq_id_dict, chain.id*string(resnum)*ins_code, "-1")),
            chain.numbering, chain.ins_codes)
        push!(chainwise_renumbering, renumbering_str)
    end

    return chainwise_renumbering
end

"""
    renumber(structure::ProteinStructure, mmcif_dict::BioStructures.MMCIFDict)

Return residue numbers from "_atom_site.label_seq_ids".

The `ProteinStructure` will automatically add a `renumbering` property if
an MMCIFDict is passed (default if the file is an MMCIF).
"""
function renumber(structure::ProteinStructure, mmcifdict)
    chainwise_renumbering = _renumber(structure, mmcifdict)
    for (i, (chain, renumbering)) in enumerate(zip(structure, chainwise_renumbering))
        structure[i] = addproperties(chain; renumbering)
    end
    return structure
end