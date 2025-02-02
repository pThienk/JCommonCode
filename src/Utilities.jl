# General utilities functions, such as reading files or miscellaneous operations
# on containers or data conversion.

"""
    Reads data within a file with the newline character as delimiter
    and return them as a Vector.

    Parameters: filename::String - REQUIRED
                type::Any - the type of data - Optional kwargs, Default=Float64

    Return: Vector{type} - the data Vector, where 'type' is chosen by the user
"""
@inline function readline_data(filename::String; type=Float64)::Vector
    return parse.(type, readlines(filename))
end

"""
    Performs 1D trapezoidal integration on a discreet function (as a Vector) with respect
    to some corresponding discreet domain (also as a Vector). The resulting Vector contains
    accumulated integration result for each domain step.
        
    Parameters: x::Vector{<:Number} - the discreet domain - REQUIRED
                f::Vector{<:Number} - the discreet function - REQUIRED

    Return: CR::Vector{Number} - cumulative integration result for each domain step
"""
@inline function cumulative_trapz_int(x::Vector{<:Real}, f::Vector{<:Real})::Vector{Real}
    CR::Vector{Real} = []

    min_indexes = min(length(x), length(f)) == length(x) ? eachindex(x) : eachindex(f)

    for i ∈ min_indexes
        push!(CR, trapz((x[1:i],), f[1:i]))
    end

    return CR
end

"""
    Finds the sliding medians of a 1D Vector according to window_size. The returned Vector
    is always of the same length as the original. The algorithm uses the nearest valid value
    for edge conditions; this is equivalent to SciPy's median_filter() with mode='nearest'.
    Currently, only this condition is supported, might be extended in the future!
    
    Parameters: samp::Vector{<:Real} - the vector of samples - REQUIRED
                window_size::Int - window size for taking the median - REQUIRED

    Return: filtered::Vector{Real} - the resulting Vector, always of the same length as samp
"""
function sliding_median(samp::Vector{<:Real}, window_size::Int)::Vector{Real}
    
end