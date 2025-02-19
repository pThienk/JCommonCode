module JCommonCode

# Imports all needed libraries
using Statistics, StatsBase
using Trapz
using Plots
using Printf

# Includes and exports all utilities functions
export readline_data
export cumulative_trapz_int
export sliding_median
export unique_count
export pow_linemaker
export scaling_plot
include("Utilities.jl")

# Includes and exports user-side slips functions
export get_slips
export get_slips_vel
include("SlipAnalysis.jl")

# Includes and exports all CCDF/CDF related functions
export ccdf
export cdf
include("CCDF.jl")

# Includes and exports all power laws MLE related functions
export mle
export ks_statistics
export scan_mle
export scan_ks
export mle_ks
include("PowerMLE.jl")





end # module JCommonCode
