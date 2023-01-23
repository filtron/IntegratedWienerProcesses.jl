module IntegratedWienerProcesses

using LinearAlgebra, SpecialFunctions

struct IWP{T<:Real,N<:Integer}
    ndiff::N
    dim::N
    IWP{T}(ndiff::N, dim::N) where {T<:Real,N<:Integer} = new{T,N}(ndiff, dim)
end


IWP(ndiff::Integer, dim::Integer) = IWP{Float64}(ndiff, dim)


"""
    ndiff(M::IWP)

Computes the number of times M is differentiable.
"""
ndiff(M::IWP) = M.ndiff

"""
    dim(M::IWP)

Computes the number of dimensions of M.
"""
dim(M::IWP) = M.dim

"""
    IWP{U}(M::IWP{T,N})

Computes an IWP of element type U from the IWP M of element type T
"""
IWP{U}(M::IWP{T,N}) where {U,T,N} = IWP{U}(ndiff(M), dim(M))



"""
    ssparams_reverse(M::IWP{T})

Computes the matrices A, B, and C associated with the state-space realization of M

dx = A*x*dt + B*dw
y = C*x

where the state vector x is ordered from highest Taylor coefficient to the lowest.
"""
function ssparams_reverse(M::IWP{T}) where {T}
    A, B, C = _iwp_ssparams_reverse_1d(T, ndiff(M))
    Id = I(dim(M))
    return kron(A, Id), kron(B, Id), kron(C, Id)
end

"""
    state2diff_matrix_reverse(M::IWP{T}, m::Integer)

Computes the matrix Cm such that

D^m y = Cm*x,

where x is the state vector of the reverse ordered state-space realization.
"""
function state2diff_matrix_reverse(M::IWP{T}, m::Integer) where {T}
    nstates = ndiff(M) + 1
    C = zeros(T, 1, nstates)
    C[1, nstates-m] = one(T)
    return kron(C, one(T) * I(dim(M)))
end

"""
    transition_parameters_cholf_reverse(M::IWP{T}, dt::T)

Computes the transition parameters A, L of the IWP M for the interval dt.
Here  L is the lower triangular Cholesky factor of the process noise,
"""
function transition_parameters_cholf_reverse(M::IWP{T}, dt::T) where {T<:Real}
    Id = one(T)I(dim(M))

    Φ1 = _transition_matrix_1d(ndiff(M), dt)
    L1 = _rgram_cholf_reverse_1d(ndiff(M), dt)

    return kron(Φ1, Id), kron(L1, Id)
end


"""
    _iwp_ssparams_reverse_1d(::Type{T}, ndiff::Integer)

Computes the matrices A, B, and C associated with the state-space realization of a one dimensiona IWP,
where the state vector x is ordered from highest Taylor coefficient to the lowest.
"""
function _iwp_ssparams_reverse_1d(::Type{T}, ndiff::Integer) where {T<:Number}
    A = diagm(-1 => ones(T, ndiff))
    B = vcat(one(T), zeros(T, ndiff, 1))
    C = hcat(zeros(T, 1, ndiff), one(T))

    return A, B, C
end

"""
    _transition_matrix_1d(ndiff::Integer, dt::T)

Computes the transition matrix A of an ndiff-fold integrator over the interval [0, dt].
"""
function _transition_matrix_1d(ndiff::Integer, dt::T) where {T<:Real}
    irange = collect(0:ndiff)
    #g = @. δt^irange / factorial(irange)
    g = @. exp(irange * log(dt) - logfactorial(irange))
    Φ = zeros(T, ndiff + one(ndiff), ndiff + one(ndiff))
    @simd ivdep for col = 0:ndiff
        @simd ivdep for row = col:ndiff
            @inbounds Φ[row+1, col+1] = g[row-col+1]
        end
    end

    return Φ
end

"""
    _rgram_cholf_reverse_1d(ndiff::Integer, dt::T)

Computes the lower triangular Cholesky factor reachability Grammian of ndiff-fold integrator over the interval [0, dt].
This is the same as the IWP process noise for a transtion over [t, t+dt].
"""
function _rgram_cholf_reverse_1d(ndiff::Integer, dt::T) where {T<:Real}
    L = zeros(T, ndiff + one(ndiff), ndiff + one(ndiff))
    @simd ivdep for n = 0:ndiff
        something = sqrt(dt) * dt^n * factorial(n)
        @simd ivdep for m = 0:ndiff
            if m <= n
                @inbounds L[n+1, m+1] =
                    something * sqrt(2 * m + 1) / factorial(n - m) / factorial(n + m + 1)
            end
        end
    end
    return L
end


export IWP,
    ndiff,
    dim,
    ssparams_reverse,
    state2diff_matrix_reverse,
    transition_parameters_cholf_reverse

end
