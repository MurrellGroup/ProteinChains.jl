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
        @test PersistentProperty(value)[1:2].value == value
        @test IndexableProperty(value)[1:2].value == value[:,1:2]
        @test (; a=PersistentProperty(1), b=IndexableProperty([1])) isa ProteinChains.NamedProperties
    end

    @testset "chain.jl" begin
        chain = ProteinChain("A", get_atoms(ProteinChains.Backbone(rand(3, 3, 5))), "AMINO", collect(1:5))
        @test length(chain) == 5
    end

    @testset "structure.jl" begin
        chain = ProteinChain("A", get_atoms(ProteinChains.Backbone(rand(3, 3, 5))), "AMINO", collect(1:5))
        structure = ProteinStructure("1CHN", Atom{Float64}[], [chain, chain])
        @test structure[1] === chain
        @test structure["A"] === chain
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
                @test ProteinChains.get_auth_asym_to_taxid(mmcifdict) == Dict("Y" => "9031", "L" => "10090", "H" => "10090")
            end
        end

    end

    @testset "store" begin
        mktempdir() do dir
            filename = joinpath(dir, "store.h5")
            structures = [addproperties(pdb"1EYE", :backbone), addproperties(pdb"3HFM", :backbone, :bond_lengths)]
            ProteinChains.serialize(filename, structures)
            structures_copy = ProteinChains.deserialize(filename)
            @test structures == structures_copy

            store = ProteinStructureStore(filename)
            @test haskey(store, "1EYE.cif")
            @test issubset((:backbone,), propertynames(store["1EYE.cif"][1]))
            @test issubset((:backbone, :bond_lengths), propertynames(store["3HFM.cif"][1]))
            delete!(store, "1EYE.cif")
            @test !haskey(store, "1EYE.cif")

            close(store)
        end
    end

end
