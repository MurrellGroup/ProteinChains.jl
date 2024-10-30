using ProteinChains
using Test

@testset "ProteinChains.jl" begin

    @testset "ideal.jl" begin
        geometry = BackboneGeometry(N_Ca_length=3, Ca_C_length=3, N_Ca_C_angle=π/2)
        ideal_residue = IdealResidue{Float64}(geometry)
        @test ideal_residue == Float64[-2 1 1; -1 -1 2; 0 0 0]
    end

    @testset "atom.jl" begin

        @testset "atom name" begin
            expected = [
                ("N", "N", " N  "),
                ("FE", "FE", "FE  "),
                ("CA", "C", " CA "),
                ("HA1", "H", " HA1"),
                ("HD11", "H", "HD11")] # observed in e.g. 5F0K
            for (name, symbol, padded) in expected
                @test ProteinChains.pad_atom_name(name, symbol) == padded
                @test ProteinChains.decode_atom_name(ProteinChains.encode_atom_name(name, symbol)) == padded
            end
        end

        @testset "Atom" begin
            coords = Float64[0, 0, 0]
            atom = Atom{Float64}(0x20414320, Int8(6), coords...)
            @test convert(Atom{Float32}, atom).x isa Float32
            @test Atom(0x20414320, 6, coords...) == atom
            @test Atom("CA", "C", coords...) == atom
            @test Atom("CA", "C", coords) == atom
            @test ProteinChains.atom_name(atom) == " CA "
            @test ProteinChains.atom_number(atom) == 6
            @test ProteinChains.atom_coords(atom) == coords
            @test ProteinChains.atom_symbol(atom) == "C"
        end

    end

    @testset "properties.jl" begin
        value = rand(2, 3)
        @test IndexableProperty(value)[1:2].value == value[:,1:2]
    end

    @testset "chain.jl" begin
        chain = ProteinChain("A", get_atoms(ProteinChains.Backbone(rand(3, 3, 5))), "AMINO", "     ", collect(1:5))
        @test length(chain) == 5
        chain = addproperties(chain, taxid=9606)
        @test hasproperty(chain, :taxid)
        chain = removeproperties(chain, :taxid)
        @test !hasproperty(chain, :taxid)
    end

    @testset "structure.jl" begin
        chain = ProteinChain("A", get_atoms(ProteinChains.Backbone(rand(3, 3, 5))), "AMINO", "     ", collect(1:5))
        structure = ProteinStructure("1CHN", Atom{Float64}[], [chain, chain])
        @test structure[1] == chain
        @test structure["A"] == chain
    end

    @testset "io" begin

        @testset "read" begin
            chains_pdb = pdbentry("1ASS"; format=PDBFormat)
            @test length.(collect(chains_pdb)) == [152]
            chains_cif = pdbentry("1ASS"; format=MMCIFFormat)
            @test chains_pdb[1].sequence == chains_cif[1].sequence
        end

        @testset "write" begin
            chains = pdb"1ASS"
            new_chains = mktempdir() do temp_dir
                temp_path = joinpath(temp_dir, "temp.pdb")
                write(temp_path, chains)
                read(temp_path, ProteinStructure)
            end
            @test chains[1].sequence == new_chains[1].sequence
        end

        @testset "mmcifutils" begin
            mktempdir() do dir
                structure = pdbentry("3HFM"; dir)
                mmcifdict = MMCIFDict(joinpath(dir, structure.name))
                auth_asym_to_taxid = mapmmcif(mmcifdict,
                    "_atom_site.auth_asym_id"   => "_atom_site.label_entity_id",
                    "_entity_src_gen.entity_id" => "_entity_src_gen.pdbx_gene_src_ncbi_taxonomy_id")
                @test auth_asym_to_taxid == Dict("Y" => "9031", "L" => "10090", "H" => "10090")
            end
        end

        # See https://proteopedia.org/wiki/index.php/Unusual_sequence_numbering
        @testset "Unusual numbering" begin
            
            @testset "Arbitrary Numbering" begin
                structure = ProteinChains.renumber(pdb"1BSZ", mmcifdict"1BSZ")
                @test !allequal(chain -> chain.numbering, structure)
                @test allequal(chain -> chain.renumbering, structure)
            end

            @testset "N-Terminal Residues Missing Coordinates" begin
                structure = ProteinChains.renumber(pdb"1D66", mmcifdict"1D66")
                chain = structure["A"]
                @test chain.numbering[begin] == 8
                @test chain.renumbering[begin] == 8
            end

            @testset "Starts With Zero Or Negative Numbers" begin
                structure = ProteinChains.renumber(pdb"1BXW", mmcifdict"1BXW")
                chain = structure["A"]
                @test chain.numbering[begin] == 0
                @test chain.renumbering[begin] == 1

                structure = ProteinChains.renumber(pdb"1D5T", mmcifdict"1D5T")
                chain = structure["A"]
                @test chain.numbering[begin] == -2
                @test chain.renumbering[begin] == 1
            end

            @testset "Insertion Codes" begin
                structure = ProteinChains.renumber(pdb"1IGY", mmcifdict"1IGY")
                chain = structure["B"]
                @test count(==(82), chain.numbering) == 4
                @test allunique(chain.renumbering)
                @test Set(collect(chain.ins_codes)) == Set([' ', 'A', 'B', 'C'])
            end

            @testset "Skipping Sequence Numbers" begin
                structure = ProteinChains.renumber(pdb"2ACE", mmcifdict"2ACE")
                chain = structure["A"]
                @test any(!isone, chain.numbering[begin+1:end] - chain.numbering[begin:end-1])
                @test chain.numbering == chain.renumbering
            end

            @testset "Not Monotonic" begin
                structure = ProteinChains.renumber(pdb"1NSA", mmcifdict"1NSA")
                chain = structure["A"]
                @test issorted(chain.numbering) # BioStructures sorts automatically with `fixlists!`, so residues are in wrong order
                @test !issorted(chain.renumbering) # the renumbering is not monotonic. it would be if BioStructures didn't sort.
            end

        end

    end

    @testset "store" begin
        mktempdir() do dir
            filename = joinpath(dir, "store.h5")
            structures = [pdb"1EYE", pdb"3HFM"]
            for (i, chain) in enumerate(structures[1])
                structures[1][i] = addproperties(chain, rand3=rand(3))
            end
            for (i, chain) in enumerate(structures[2])
                structures[2][i] = addproperties(chain, rand3=rand(3), taxid=-1)
            end
            ProteinChains.serialize(filename, structures)
            structures_copy = ProteinChains.deserialize(filename)
            @test structures == structures_copy

            store = ProteinStructureStore(filename)
            @test haskey(store, "1EYE.cif")
            @test issubset((:rand3,), propertynames(store["1EYE.cif"][1]))
            @test all(chain -> issubset((:rand3, :taxid), propertynames(chain)), store["3HFM.cif"])
            delete!(store, "1EYE.cif")
            @test !haskey(store, "1EYE.cif")

            close(store)
        end
    end

end
