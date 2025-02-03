module JCommonCode

# Imports all needed libraries
using Statistics, StatsBase
using Trapz

# Includes and exports all utilities functions
export readline_data
export cumulative_trapz_int
export sliding_median
export unique_count
export pow_linemaker
include("Utilities.jl")

# Includes and exports user-side slips functions
export get_slips
export get_slips_vel
include("SlipAnalysis.jl")

# Includes and exports all CCDF/CDF related functions
export ccdf
export cdf
include("CCDF.jl")





end # module JCommonCode
