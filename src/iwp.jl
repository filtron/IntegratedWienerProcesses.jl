struct IWP{T,N<:Integer} <: AbstractIWP{T}
    ndiff::N
    dim::N
    IWP{T}(ndiff::N, dim::N) where {T<:Number,N<:Integer} = new{T,N}(ndiff, dim)
end


IWP(ndiff::Integer, dim::Integer) = IWP{Float64}(ndiff, dim)


"""
    ndiff(M::IWP)

Computes the number of times M is differentiable.
"""
ndiff(model::IWP) = model.ndiff

"""
    dim(M::IWP)

Computes the number of dimensions of M.
"""
dim(model::IWP) = model.dim

"""
    ssparams(model::IWP{T}, method::AbstractRealizationMethod) where {T}

Computes the state-space matrices A, B, C for model model using realization method method.
"""
function ssparams(model::IWP{T}, method::AbstractRealizationMethod) where {T}
    A, B, C = _iwp_ssparams_1d(T, ndiff(model), method)
    Id = I(dim(model))
    return kron(A, Id), kron(B, Id), kron(C, Id)
end

function _iwp_ssparams_1d(::Type{T}, ndiff::Integer, ::ReverseTaylor) where {T<:Number}
    A = diagm(-1 => ones(T, ndiff))
    B = vcat(one(T), zeros(T, ndiff, 1))
    C = hcat(zeros(T, 1, ndiff), one(T))
    return A, B, C
end


"""
    state2diff_matrix(model::IWP{T}, m::Integer, method::AbstractRealizationMethod) where {T}

Computes the matrix Em such that Dy = Em * x, where x is the state vector associated with realization method method.
"""
function state2diff_matrix(model::IWP{T}, m::Integer, ::ReverseTaylor) where {T}
    nstates = ndiff(model) + 1
    C = zeros(T, 1, nstates)
    C[1, nstates-m] = one(T)
    return kron(C, one(T) * I(dim(model)))
end

"""
    transition_matrix(model::IWP{T}, dt::Real, method::AbstractRealizationMethod) where {T}

Computes the transition matrix of model model, using step-size dt, and realization method method.
"""
function transition_matrix(
    model::IWP{T},
    dt::Real,
    method::AbstractRealizationMethod,
) where {T}
    V = promote_type(real(T), typeof(dt))
    return kron(_transition_matrix_1d(ndiff(model), V(dt), method), one(T)I(dim(model)))
end

function _transition_matrix_1d(ndiff::Integer, dt::T, ::ReverseTaylor) where {T<:Real}
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
    transition_cov_cholf(model::IWP{T}, dt::T, method::AbstractRealizationMethod) where {T}

Computes a left cholesky factor of the transition covariance of the model model, with step-size dt, using realization method method.
"""
function transition_cov_cholf(
    model::IWP{T},
    dt::Real,
    method::AbstractRealizationMethod,
) where {T}
    V = promote_type(real(T), typeof(dt))
    return kron(transition_cov_cholf_1d(ndiff(model), V(dt), method), one(T)I(dim(model)))
end

function transition_cov_cholf_1d(ndiff::Integer, dt::T, ::ReverseTaylor) where {T}
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

"""
    transition_cov(model::IWP{T}, dt::T, method::AbstractRealizationMethod) where {T}

Computes the transition covariance of the model model, with step-size dt, using realization method method.
"""
function transition_cov(
    model::IWP{T},
    dt::Real,
    method::AbstractRealizationMethod,
) where {T}
    V = promote_type(real(T), typeof(dt))
    return kron(transition_cov_1d(ndiff(model), V(dt), method), one(T)I(dim(model)))
end

function transition_cov_1d(ndiff::Integer, dt::T, ::ReverseTaylor) where {T}
    Q = zeros(ndiff + 1, ndiff + 1)
    @simd ivdep for n = 0:ndiff
        @simd ivdep for m = 0:ndiff
            Q[n+1, m+1] = dt^(n + m + 1) / (n + m + 1) / factorial(n) / factorial(m)
        end
    end
    return Q
end
