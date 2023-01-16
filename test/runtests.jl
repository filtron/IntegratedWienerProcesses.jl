using IntegratedWienerProcesses, LinearAlgebra
using Test

@testset "IntegratedWienerProcesses.jl" begin

    @testset "IWP" begin

        d = 2
        M = IWP(1, d)

        dt = 0.5

        A = kron([1.0 0.0; dt 1.0], 1.0I(d))
        L = kron([sqrt(dt) 0.0; sqrt(dt^3)/2 sqrt(dt^3 / 12.0)], 1.0I(d))

        @test ndiff(M) == 1
        @test state2diff_matrix_reverse(M, 0) ==
              kron(hcat(zeros(1, 1), fill(1.0, 1, 1)), 1.0I(d))
        @test all(transition_parameters_cholf_reverse(M, dt) .≈ (A, L))
    end

end
