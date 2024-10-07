using ProteinChains
using Test

@testset "ProteinChains.jl" begin

    @testset "ProteinChain" begin
        chain = ProteinChain("A", get_atoms(ProteinChains.Backbone(rand(3, 3, 5))), "AMINO", collect(1:5))
        @test length(chain) == 5
    end

    @testset "ProteinStructure" begin
        chain = ProteinChain("A", get_atoms(ProteinChains.Backbone(rand(3, 3, 5))), "AMINO", collect(1:5))
        structure = ProteinStructure("1CHN", Atom{Float64}[], [chain])
        @test structure[1] === chain
        @test structure["A"] === chain
    end

    @testset "read/write" begin

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

    end

    @testset "serialization" begin
        mktempdir() do dir
            filename = "dataset.h5"
            path = joinpath(dir, filename)
            dataset1 = ProteinDataset([annotate(pdb"1EYE", :backbone)])
            ProteinChains.serialize(path, dataset1)
            dataset2 = ProteinChains.deserialize(path)
            @test dataset1 == dataset2
            @test hasproperty(first(values(dataset1))["A"], :backbone)
        end
    end

end
