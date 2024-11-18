using ProteinChains
using Test

@testset "ProteinChains.jl" begin

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

    end

    @testset "store" begin
        mktempdir() do dir
            filename = joinpath(dir, "store.h5")
            structures = [pdb"1EYE", pdb"3HFM"]
            for (i, chain) in enumerate(structures[1])
                addproperties!(chain, rand3=rand(3))
            end
            for (i, chain) in enumerate(structures[2])
                addproperties!(chain, rand3=rand(3), taxid=-1, dict=Dict("a"=>1, "b"=>2))
            end
            ProteinChains.serialize(filename, structures)
            structures_copy = ProteinChains.deserialize(filename)
            @test structures == structures_copy

            store = ProteinStructureStore(filename)
            @test haskey(store, "1EYE.cif")
            @test issubset((:rand3,), propertynames(store["1EYE.cif"][1]))
            @test all(chain -> issubset((:rand3, :taxid, :dict), propertynames(chain)), store["3HFM.cif"])
            @test store["3HFM.cif"][1].dict == Dict("a"=>1, "b"=>2)
            delete!(store, "1EYE.cif")
            @test !haskey(store, "1EYE.cif")

            close(store)
        end
    end

end
