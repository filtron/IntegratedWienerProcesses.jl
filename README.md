# IntegratedWienerProcesses

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://filtron.github.io/IntegratedWienerProcesses.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://filtron.github.io/IntegratedWienerProcesses.jl/dev/)
[![Build Status](https://github.com/filtron/IntegratedWienerProcesses.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/filtron/IntegratedWienerProcesses.jl/actions/workflows/CI.yml?query=branch%3Amain)

A package implementing tools for working with integrated Wiener processes. 
A $\nu$-times integrated wiener process on the interval $[0, T]$ is given by 

$$
y(t) = \sum_{k=0}^\nu y^{(k)}(0) \frac{t^k}{k!} + \int_0^t \frac{(t-\tau)^\nu}{\nu!} \operatorname{d} w(t), 
$$

where $y^{(k)}$ is the $k$:th derivative of $y$ and $w$ is a standard Wiener process. 
Integrated Wiener processes admit state-space realizations according to 

$$
\begin{align}
\operatorname{d} x(t) &= A x(t) \operatorname{d} t + B \operatorname{d} w(t), \\ 
y(t) &= C x(t). 
\end{align}
$$
