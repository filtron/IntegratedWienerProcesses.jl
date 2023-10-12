using IntegratedWienerProcesses, LinearAlgebra
using Test

@testset "IntegratedWienerProcesses.jl" begin

    dt = 0.5
    method = ReverseTaylor()

    @testset "IWP" begin

        d = 2
        model = IWP(1, d)

        A = kron([1.0 0.0; dt 1.0], 1.0I(d))
        L = kron([sqrt(dt) 0.0; sqrt(dt^3)/2 sqrt(dt^3 / 12.0)], 1.0I(d))

        @test ndiff(model) == 1
        @test state2diff_matrix(model, 0, method) ==
              kron(hcat(zeros(1, 1), fill(1.0, 1, 1)), 1.0I(d))
        @test transition_matrix(model, dt, method) ≈ A
        @test transition_cov_cholf(model, dt, method) ≈ L
    end

    @testset "transition_cov_cholf" begin

        for i = 0:9
            m = IWP(i, 1)
            L = transition_cov_cholf(m, dt, method)
            @test L * L' ≈ transition_cov(m, dt, method)
        end

    end

    @testset "preconditioning" begin 
        for ndiff = 0:9
            L = IntegratedWienerProcesses.transition_cov_cholf_1d(ndiff, dt, ReverseTaylor())
            precond, Lbreve = IntegratedWienerProcesses.transition_cov_cholf_precond_1d(ndiff, dt, ReverseTaylor())
            @test L ≈ Diagonal(precond) * Lbreve
        end
    end


end
