# General utilities functions, such as reading files or miscellaneous operations
# on containers or data conversion.

"""
    Read the data within a file with the newline character as delimiter.

    return Vector{T<:Number} where T is chosen by the user; Default=Float64
"""
function readline_data(filename::String; type=Float64)::Vector{<:Number}
    return parse.(type, readlines(filename))
end