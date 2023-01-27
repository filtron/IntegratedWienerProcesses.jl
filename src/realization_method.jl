

abstract type AbstractRealizationMethod end

"""
    ReverseTaylor

Type for representing state-space realizations in terms of Taylor coefficients in reverse order.
That is the state-vector is of the form:

    x = [D^n y, D^(n-1)y, ..., Dy, y]
"""
struct ReverseTaylor <: AbstractRealizationMethod end
