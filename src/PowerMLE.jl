# Implementation of the Maximum Likelihood Estimation of power laws and Kolmogorov-Smirnov statistics
# described in 
# Clauset, A., Shalizi, C. R., & Newman, M. E. J. (2009).
# Power-Law Distributions in Empirical Data. SIAM Review, 51(4),
# 661–703. http://dx.doi.org/10.1137/070710111, https://arxiv.org/abs/0706.1062

# Parts modified from https://github.com/jlapeyre/MaximumLikelihoodPower.jl by John Lapeyre to be compatible
# with newer Julia versions.

"""
    MLEScan{T <: Real}

    Record best estimate of alpha and associated parameters.
"""
mutable struct MLEScan
    alpha::Real
    stderr::Real
    minKS::Real
    xmin::Real
    imin::Int
    npts::Int
    nptsall::Int
    ntrials::Int
end

"""
    MLEKS{T}

    Container for storing results of MLE estimate and
    Kolmogorov-Smirnov statistics of the exponent of a power law.
"""
struct MLEKS{T}
    alpha::T
    stderr::T
    KS::T
end

function MLEScan(T)
    z = zero(T)
    return MLEScan(z, z, Inf, z, 0, 0, 0, 0)
end

"""
    Overloads the standard show() function to prints out the state of MLEScan
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
    Returns the maximum likelihood estimate and standard error of the exponent of a power law
    applied to a Vector.

    Parameter: data::AbstractVector{<:Real} - Vector of data points - REQUIRED

    Returns: Tuple of the MLE and SE in the form (ahat::Real, stderr::Real)
"""
function mle(data::AbstractVector{<:Real})::Tuple{Real, Real}

    @assert (!isempty(data)) "The data Vector cannot be empty!"

    data = sort(data; alg=QuickSort)
    
    xmin::Real = data[1]
    acc::Real = 0
    xlast::Real = Inf
    ncount::Int = 0
    @inbounds for value ∈ data
        xlast == value && continue
        xlast = value
        ncount += 1
        acc += log(value / xmin)
    end

    ahat::Real = 1 + ncount / acc
    stderr::Real = (ahat - 1) / sqrt(ncount)
    return (ahat, stderr)
end

"""
    Returns the Kolmogorov-Smirnov statistics (max distance) comparing data to a power law with power alpha.
    
    Parameters: data::AbstractVector{<:Real} - the Vector of data points - REQUIRED
                alpha::Real - the power to compare - REQUIRED

    Return: max_distance::Real - the Kolmogorov-Smirnov statistics max distance
"""
function ks_statistics(data::AbstractVector{<:Real}, alpha::Real)
    
    @assert (!isempty(data)) "The data Vector cannot be empty!"

    data = sort(data; alg=QuickSort)
    unique!(data)

    num::Int = length(data)
    xmin::Real = data[1]
    max_distance::Real = 0
    @inbounds for i ∈ 0:num-1
        pl::Real = 1 - (xmin / data[i+1])^alpha
        distance::Real = abs(pl - i / num)

        if distance > max_distance
            max_distance = distance
        end
    end

    return max_distance
end

"""
    Computes the Kolmogorov Smirnov statistics for several values of α in the iterator powers.

    Parameters: data::AbstractVector{<:Real} - the Vector of data points - REQUIRED
                powers::AbstractVector{<:Real} - the iterator of powers to test - REQUIRED

    Returns: the value of α that minimizes the KS statistic and the two neighboring values.
"""
function scan_ks(data::AbstractVector{<:Real}, powers::AbstractVector{<:Real})

    @assert (!isempty(data)) "The data Vector cannot be empty!"

    ks::Vector{Real} = [ks_statistics(data, p) for p in powers]
    i = argmin(ks)
    return powers[(i-1):(i+1)] |> collect
end

"""
    Returns the MLE and SE of the exponent of a power law
    applied to the sorted Vector data. Also return the Kolmogorov-Smirnov statistics.
    
    Parameter: data::AbstractVector{<:Real} - the Vector of data points - REQUIRED
    
    Return: results are returned in an instance of type MLEKS.
"""
function mle_ks(data::AbstractVector{<:Real})

    @assert (!isempty(data)) "The data Vector cannot be empty!"

    alpha, stderr = mle(data)
    KSstat = ks_statistics(data, alpha)
    return MLEKS(alpha, stderr, KSstat)
end

"""
    Copies and mutate MLEScan object
"""
function copy_mslescan!(mlescan::MLEScan, mle::MLEKS, data::AbstractVector{<:Real}, i::Integer)

    mlescan.minKS = mle.KS
    mlescan.alpha = mle.alpha
    mlescan.stderr = mle.stderr
    mlescan.imin = i
    mlescan.npts = length(data)
    mlescan.xmin = data[1]
end

"""
    Compares the results of MLE estimation to record results in mlescan and update mlescan.
"""
function compare_scan(mlescan::MLEScan, mle::MLEKS, data::AbstractVector{<:Real}, i::Integer)

    if mle.KS < mlescan.minKS
        copy_mslescan!(mlescan, mle, data, i)
    end
    mlescan.ntrials += 1

end

"""
    Performs MLE approximately ntrials times on data, increasing xmin. Stop trials
    if the standard error of the estimate alpha is greater than stderrcutoff.
    If useKS is true, then the application of MLE giving the smallest KS statistics is
    returned. Returns an object containing statistics of the scan.

    scan_mle is intended to analayze the power-law behavior of the tail of data.

    Parameters: data::AbstractVector{<:Real} - the Vector of data points - REQUIRED
                ntrials::Int - number of trials to perform - Optional; Default=100
                stderrcutoff::Real - maximum threshold SE to stop trials - Optional; Default=0.1
                useKS::Bool - whether to use the smallest MLE - Optional; Default=false

    Returns: mlescan::MLEScan - the object containing scans' results
"""
function scan_mle(data::AbstractVector{<:Real}; ntrials::Int=100, stderrcutoff::Real=0.1, useKS::Bool=false)

    @assert (!isempty(data)) "The data Vector cannot be empty!"

    data = sort(data; alg=QuickSort)

    skip::Int = round(Int, length(data) / ntrials)
    if skip < 1
        skip = 1 
    end

    return _scan_mle(data, 1:skip:length(data), stderrcutoff, useKS)
end

"""
    Backend of scan_mle
"""
function _scan_mle(data, range::AbstractVector{<:Integer}, stderrcutoff, useKS)

    mlescan = MLEScan(eltype(data))
    mlescan.nptsall = length(data)
    lastind::Int = 0
    @inbounds for i in range
        ndata = @view data[i:end]
        mleks = mle_ks(ndata)
        lastind = i
        if mleks.stderr > stderrcutoff || i == last(range)
            if !useKS
                copy_mslescan!(mlescan, mleks, ndata, i)
                mlescan.ntrials = i
            end
            break
        end

        if useKS
            compare_scan(mlescan, mleks, ndata, i)
        end  # do we want ndata or data here ?

    end

    return mlescan
end