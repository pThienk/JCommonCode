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
    is always of the same length and type as the original. The algorithm uses nearest valid values
    for edge conditions - this is equivalent to SciPy's median_filter() with mode='nearest'.
    Currently, only this condition is supported; might be extended in the future!
    
    Parameters: samp::Vector{<:Real} - the vector of samples - REQUIRED
                window_size::Int - window size for taking the median, should be odd - REQUIRED

    Return: filtered::Vector{eltype(samp)} - the resulting Vector, always of the same length and type as samp
"""
@inline function sliding_median(samp::Vector{<:Real}, window_size::Int)
    @assert (window_size <= length(samp)) "Window size should not be larger than sample size!"
    @assert isodd(window_size) "Window size should be odd for maximum efficiency!"

    filtered::Vector{eltype(samp)} = []
    half_win_size::Int = (window_size ÷ 2) + 1

    for i ∈ eachindex(samp)

        window_vec::Vector{Real} = []
        if i < half_win_size
            edge = [samp[1] for _ in 1:(half_win_size - i)]
            prepend!(window_vec, edge)
            append!(window_vec, samp[1:i])
        elseif i + half_win_size - 1 > lastindex(samp)
            edge = [samp[end] for _ in 1:(i + half_win_size - 1 - lastindex(samp))]
            prepend!(window_vec, samp[i:end])
            append!(window_vec, edge)
        else
            append!(window_vec, samp[(i-half_win_size+1):(i+half_win_size-1)])
        end

        push!(filtered, median(window_vec))
    end

    return filtered
end