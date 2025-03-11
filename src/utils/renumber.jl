_nothing_fallback(x, y) = x
_nothing_fallback(::Nothing, y) = y
tryparse_int32(s::String) = _nothing_fallback(tryparse(Int32, s), Int32(-1))

function renumber(chain::ProteinChain, mmcif_dict::BioStructures.MMCIFDict)
    id = split(chain.id, '-')[1]

    label_seq_ids = mmcif_dict["_atom_site.label_seq_id"]
    pdbx_PDB_ins_codes = mmcif_dict["_atom_site.pdbx_PDB_ins_code"]
    auth_seq_ids = mmcif_dict["_atom_site.auth_seq_id"]
    auth_asym_ids = mmcif_dict["_atom_site.auth_asym_id"]

    label_seq_id_dict = Dict{Tuple{Int32,Char},String}()

    current_residue = ("", ' ')
    for (i, auth_asym_id) in enumerate(auth_asym_ids)
        auth_asym_id == id || continue
        ins_code = pdbx_PDB_ins_codes[i] == "?" ? " " : pdbx_PDB_ins_codes[i]

        residue = (auth_seq_ids[i], first(ins_code))
        if current_residue == residue
            continue
        else
            current_residue = residue
        end

        renum = tryparse_int32(auth_seq_ids[i])

        key = (renum, first(ins_code))
        label_seq_id_dict[key] = label_seq_ids[i]
    end

    renumbering = map(
        (num, ins_code) -> tryparse_int32(get(label_seq_id_dict, (num, ins_code), "-1")),
        chain.numbering, chain.ins_codes)

    return renumbering
end

function renumber!(chain::ProteinChain, mmcif_dict::BioStructures.MMCIFDict)
    chain.numbering = Int.(renumber(chain, mmcif_dict))
    return chain
end
