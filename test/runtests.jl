using ProteinChains
using Test

@testset "ProteinChains.jl" begin

    @testset "ProteinChain" begin
        chain = ProteinChain("A", "AMINO", rand(3, 3, 5), collect(1:5), [ProteinChains.Atom{Float64}[] for _ in 1:5])
        @test length(chain) == 5
        @test summary(chain) == "5-residue ProteinChain{Float64} (A)"

        i = 5:-1:1
        @test chain[i].backbone == chain.backbone[:,:,i]
    end

    @testset "ProteinStructure" begin
        chain = ProteinChain("A", "AMINO", rand(3, 3, 5), collect(1:5), [ProteinChains.Atom{Float64}[] for _ in 1:5])
        structure = ProteinStructure("1CHN", [chain])
        @test structure[1] === chain
        @test structure["A"] === chain
    end

    @testset "AnnotatedProteinChain" begin
        chain = ProteinChain("A", "AMINO", rand(3, 3, 5), collect(1:5), [ProteinChains.Atom{Float64}[] for _ in 1:5])
        annotated_chain = annotate(chain; modelnum=1)
        @test annotated_chain isa AnnotatedProteinChain
        @test annotated_chain.modelnum == 1
        annotate_indexable!(annotated_chain; secondary_structure=assign_secondary_structure(annotated_chain))
        @test length(annotated_chain.secondary_structure) == 5
        @test length(annotated_chain[1:3].secondary_structure) == 3
    end

    @testset "read/write" begin

        @testset "read" begin
            chains_pdb = pdbentry("1ASS"; format=PDBFormat)
            @test length.(collect(chains_pdb)) == [152]
            chains_cif = pdbentry("1ASS"; format=MMCIFFormat)
            @test chains_pdb[1].backbone == chains_cif[1].backbone
        end

        @testset "write" begin
            chains = pdb"1ASS"
            new_chains = mktempdir() do temp_dir
                temp_path = joinpath(temp_dir, "temp.pdb")
                write(temp_path, chains)
                read(temp_path, ProteinStructure)
            end
            @test chains[1].backbone == new_chains[1].backbone
        end

    end

end
