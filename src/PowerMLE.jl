# Implementation of the Maximum Likelihood Estimation of power laws and Kolmogorov-Smirnov statistics
# described in 
# Clauset, A., Shalizi, C. R., & Newman, M. E. J. (2009).
# Power-Law Distributions in Empirical Data. SIAM Review, 51(4),
# 661–703. http://dx.doi.org/10.1137/070710111, https://arxiv.org/abs/0706.1062

"""
    MLEScan{T <: Real}

    Record best estimate of alpha and associated parameters.
"""
mutable struct MLEScan{T <: Real}
    alpha::T
    stderr::T
    minKS::T
    xmin::T
    imin::Int
    npts::Int
    nptsall::Int
    ntrials::Int
end

"""
    MLEKS{T}

    Container for storing results of MLE estimate and
    Kolmogorov-Smirnov statistic of the exponent of a power law.
"""
struct MLEKS{T}
    alpha::T
    stderr::T
    KS::T
end

"""
    Overload the standard show() function to prints out the state of MLEScan
"""
function Base.show(io::IO, s::MLEScan)
    @printf(io, "alpha   = %.8f\n" , s.alpha)
    @printf(io, "stderr  = %.8f\n" , s.stderr)
    println(io, "minKS   = ", s.minKS)
    println(io, "xmin    = ", s.xmin)
    println(io, "imin    = ", s.imin)
    println(io, "npts    = ", s.npts)
    println(io, "nptsall = ", s.nptsall)
    @printf(io, "pct pts = %.3f\n", (s.npts / s.nptsall))
    println(io, "ntrials = ", s.ntrials)
    return nothing
end

"""
    
"""
function mle(data::Vector{<:Real})::Tuple{Real, Real}

    data = sort(data; alg=QuickSort)
    
    xmin::Real = data[1]
    acc::Real = 0
    xlast::Real = Inf
    ncount::Int = 0
    for value ∈ data 
        xlast == value && continue
        xlast = value
        ncount += 1
        acc += log(value / xmin)
    end

    ahat::Real = 1 + ncount / acc
    stderr::Reat = (ahat - 1) / sqrt(ncount)
    return (ahat, stderr)
end