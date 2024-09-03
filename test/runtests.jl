using ProteinChains
using Test

@testset "ProteinChains.jl" begin

    @testset "ProteinChain" begin
        chain = ProteinChain("A", "AMINO", rand(3, 3, 5), [ProteinChains.Atom{Float64}[] for _ in 1:5]; modelnum=1)
        @test countresidues(chain) == 5
        @test summary(chain) == "5-residue ProteinChain \"A\""
    end

    @testset "ProteinStructure" begin
        chain = ProteinChain("A", "AMINO", rand(3, 3, 5), [ProteinChains.Atom{Float64}[] for _ in 1:5])
        structure = ProteinStructure("1CHN", [chain])
        @test structure[1] === chain
        @test structure["A"] === chain
    end

    @testset "read/write" begin

        @testset "read" begin
            chains_pdb = pdbentry("1ASS"; format=PDBFormat)
            @test countresidues.(collect(chains_pdb)) == [152]
            chains_cif = pdbentry("1ASS"; format=MMCIFFormat)
            @test chains_pdb[1].backbone == chains_cif[1].backbone
        end

        @testset "write" begin
            chains = pdb"1ASS"
            new_chains = mktempdir() do temp_dir
                temp_path = joinpath(temp_dir, "temp.pdb")
                writechains(temp_path, chains)
                readchains(temp_path)
            end
            @test chains[1].backbone == new_chains[1].backbone
        end

    end

    @testset "secondary structure" begin
        ss_composition(chain::ProteinChain) = [count(==(ss), chain.secondary_structure) for ss in ['-', 'H', 'E']]
        structure = assign_secondary_structure!(pdb"1ZAK")
        @test ss_composition.(collect(structure)) == [[72, 116, 32], [72, 116, 32]]
    end

end
