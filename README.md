# IntegratedWienerProcesses

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://filtron.github.io/IntegratedWienerProcesses.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://filtron.github.io/IntegratedWienerProcesses.jl/dev/)
[![Build Status](https://github.com/filtron/IntegratedWienerProcesses.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/filtron/IntegratedWienerProcesses.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/filtron/IntegratedWienerProcesses.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/filtron/IntegratedWienerProcesses.jl)


A package implementing tools for working with integrated Wiener processes. 
A $\nu$-times integrated wiener process on the interval $[0, T]$ is given by 

$$
y(t) = \sum_{k=0}^\nu y^{(k)}(0) \frac{t^k}{k!} + \int_0^t \frac{(t-\tau)^\nu}{\nu!} \mathrm{d} w(t), 
$$

where $y^{(k)}$ is the $k$:th derivative of $y \in \mathbb{R}^d$ and $w$ is a standard Wiener process on $\mathbb{R}^d$. 
Integrated Wiener processes admit state-space realizations according to 

$$
\begin{align}
\operatorname{d} x(t) &= A x(t) \mathrm{d} t + B \mathrm{d} w(t), \\ 
y(t) &= C x(t). 
\end{align}
$$

The state vector is Gauss-Markov with transition distribution given by

$$
x(t+s) \mid x(t) \sim  \mathcal{N}( e^{A s} x(t), Q(s) ), 
$$

where $e^{As}$ is the transition matrix and $Q(s)$ is the process noise covariance. 
Since $Q(s)$ is positive definite it admits a Cholesky factorization $Q(s) = L(s) L^*(s)$. 
This pacakge implements functionality to represent integrated Wiener processes and work with their state-space realization. 


## Functionality 

The package defines the following type for representing $\nu$-times integrated Wiener processes. 

```julia
IWP{T,N} where {T<:Real,N<:Integer}
```
An instance may be constructed according to

```julia
IWP{T}(ndiff::Integer, dim::Integer) # ndiff times differentiable IWP of dimension dim and element type T
IWP(ndiff::Integer, dim::Integer)    # same as IWP{Float64}(ndiff, dim)
```

Additionally, the following functions are implemented. 

```julia
ndiff(M::IWP) # number of times M can be differentiated 
dim(M::IWP)   # dimension of M 
ssparams_reverse(M::IWP{T}) # computes the matrices A, B, C when the state vector is the Taylor coefficients in reverse order
state2diff_matrix_reverse(M::IWP{T}, m::Inter) # computes a matrix such that when mutiplied with the state vector the mth derivative is obtained
transition_parameters_cholf_reverse(M::IWP{T}, dt::T) # computes the transition matrix and the Cholesky factor of the process nosie covariance
```

