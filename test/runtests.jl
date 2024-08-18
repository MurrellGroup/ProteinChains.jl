using ProteinChains
using Test

@testset "ProteinChains.jl" begin

    @testset "ProteinChain" begin
        chain = ProteinChain("A", "AMINO", rand(3, 3, 5); modelnum=1)
        @test countresidues(chain) == 5
        @test summary(chain) == "5-residue ProteinChain \"A\" with 1 properties"
    end

    @testset "ProteinStructure" begin
        chain = ProteinChain("A", "AMINO", rand(3, 3, 5))
        structure = ProteinStructure("1CHN", [chain])
        @test structure[1] === chain
        @test structure["A"] === chain
        
    end

    @testset "PDB" begin

        @testset "read" begin
            structure = pdb"1ASS"
            @test countresidues.(chains_pdb) == [152]
            chains_cif = readchains("data/1ASS.cif", MMCIFFormat)
            @test chains_pdb[1].backbone.coords == chains_cif[1].backbone.coords
        end

        @testset "write" begin
            chains = pdb"1ASS"
            new_chains = mktempdir() do temp_dir
                temp_path = joinpath(temp_dir, "temp.pdb")
                writepdb(temp_path, chains)
                readpdb(temp_path)
            end
            @test chains[1].backbone.coords == new_chains[1].backbone.coords
        end

    end

end
