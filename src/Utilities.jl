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
    Read data from an entire directory of files with the specified prefix and suffix,
    else attempt to read all files in the directory. Data is returned as a vector of 
    tuples with elements of the form (filename::String, Vector{type}). 
    
    WARNING!: There is no check for the validity of the files' format. Please make sure that all files you
    are attempting to read have the valid format!
   
    Parameter:  dir::String - path to the directory to read - REQUIRED
                pref::String - prefix identifier of the files to read - Optional; Default=empty
                suff::String - suffix identifier of the files to read, INCLUDING EXTENSIONS - Optional; Default=empty
                type - type of data to be read, passes to readline_data - Optional; Default=Float64

    Return: data_bundle::Vector{Tuple} - a Vector of (filename, data::Vector{type}) Tuple, where data corresponds
            to the filename it was read from.

"""
@inline function readline_data_bundle(dir::String; pref::String="", suff::String="", type=Float64)::Vector{Tuple}

    files::Vector{String} = readdir(dir)
    dir = dir[1:end-length((split(dir, "/")[end]))]

    if !isempty(pref) && !isempty(suff)
        files = [f for f in files if (f[1:lastindex(pref)] == pref) && (f[end-lastindex(suff)+1:end] == suff)]
    elseif !isempty(pref)
        files = [f for f in files if (f[1:lastindex(pref)] == pref)]
    elseif !isempty(suff)
        files = [f for f in files if (f[end-lastindex(suff)+1:end] == suff)]
    else

    end

    data_bundle::Vector{Tuple} = []
    @inbounds for file ∈ files
        push!(data_bundle, (file, readline_data(dir * file; type=Float64)))
    end
    
    return data_bundle
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

    @inbounds for i ∈ min_indexes
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

    @inbounds for i ∈ eachindex(samp)

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

"""
    Takes a Vector and and modified it by sorting and removing repetitive elements.
    The function returns the count of each element's repetitions.
    
    WARNING!: This function will modify the original by sorting it and removing repetitions!

    Parameter: vec::Vector{<:Real} - the operating Vector - REQUIRED

    Return: Vector{Int} - a Vector of the counts of each element in the same order as 
            the Vector elements
"""
@inline function unique_count!(vec::Vector{<:Real})::Vector{Int}
    
    sort!(vec; alg=QuickSort)
    counts::Vector{Int} = []

    count::Int = 1
    current = nothing
    next = nothing
    @inbounds for i in eachindex(vec) 
        
        current = vec[i]
        next = i < lastindex(vec) ? vec[i+1] : nothing

        if current == next
            count += 1
        else
            push!(counts, count)
            count = 1
        end
    end

    unique!(vec)

    return counts
end

"""
    Returns X and Y arrays of a power-law line.

    Parameters
    ----------
    slope: (Float; REQUIRED)
        Power law PDF slope
    intercept: Tuple
        Intercept of the line
        Formatted as (x, y)
    xmin: (Float; REQUIRED)
        Minimum x-value the line will appear over
    xmax: (Float; REQUIRED)
        Maximum x-value the line will appear over
    ppd: (Int; optional)
        Number of log-spaced points per decade to evaluate the line at

    Returns
    -------
    [0] x_vals: (Array) X values of the line
    [1] y_vals: (Array) Y values of the line
"""
@inline function pow_linemaker(slope::Real, intercept::Tuple{<:Real,<:Real}, xmin::Real, xmax::Real; ppd::Real=40)::Tuple
    
    log_x_intercept, log_y_intercept = log10.(intercept)
    log_xmin = log10(xmin)
    log_xmax = log10(xmax)

    log_b = log_y_intercept - slope * log_x_intercept

    #num = round(ppd * (log_xmax - log_xmin))
    #spacing = (log_xmax - log_xmin) / num

    x_values = logrange(xmin, xmax, round(Int, ppd * (log_xmax - log_xmin))) # [10^expo for expo in log_xmin:spacing:log_xmax]
    y_values = (10^log_b) .* (x_values .^ slope)

    return (x_values, y_values)
end

"""
    Makes a simple CCDF scatter plot of the data provided. Additionally, the \'comparable\'
    parameter can be provided to plot a comparison curve alongside data. If \'save\' is provided,
    saves the plot figure at the specified path. The default scale is :log10, but can be chosen by
    the user. Note: this functionality is intended for data visualization only, not for
    publication-quality plots!

    Parameters: list_x::Vector{<:Real}, list_y::Vector{<:Real} - x and y Vectors of coordinates of the data - REQUIRED
                scale::Symbol - scale of x and y axises - Optional; Default = :log10
                save::String - path to save plot figure, does not save if empty - Optional; Default = empty
                comparable::Tuple - must be provided in the form (lx::Vector{<:Real}, ly::Vector{<:Real}, name::String),
                                    plots a comparison curve alongside data - Optional; Default = empty

    Return: the Plots object \'p\'
"""
@inline function scaling_plot(list_x::Vector{<:Real}, list_y::Vector{<:Real}; scale::Symbol=:log10, save::String="", comparable::Tuple=())

    @assert (length(list_x) == length(list_y)) "Incompatable Lists: list_x and list_y must have the same length!"

    if scale == :log10
        list_x = [x for x in list_x if (x > 0)]
        list_y = [y for y in list_y if (y > 0)]
    end

    p = scatter(list_x, list_y; label="data", color=:blue, xlabel="Sizes (s)", ylabel="P (S > s)", title="Plot of CCDF (Event Probs vs Event Sizes)",
     xscale=scale, yscale=scale, legend=:outertop, legendcolumns=(isempty(comparable) ? 1 : 2))

    if !isempty(comparable)
        lx, ly, name = comparable

        @assert (length(lx) == length(ly)) "Incompatable Lists: lx and ly in comparable must have the same length!"
        plot!(p, lx, ly; label=name, color=:red, lw=3)
    end

    if !isempty(save)
        savefig(p, save)
    end

    gui()
    return p
end