module IntegratedWienerProcesses

using LinearAlgebra, SpecialFunctions

include("realization_method.jl")
export AbstractRealizationMethod, ReverseTaylor

abstract type AbstractIWP{T<:Number} end

include("iwp.jl")
export AbstractIWP,
    IWP,
    ndiff,
    dim,
    ssparams,
    state2diff_matrix,
    transition_matrix,
    transition_cov_cholf,
    transition_cov

end
