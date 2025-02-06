include("../src/JCommonCode.jl")

using .JCommonCode
using Plots

stress = readline_data("test/test_data/test_data_2.txt")
dstress = diff(stress)

v, t, s, d, b, e = get_slips(; vel=dstress, drops=true, threshold=-1)

x, y = ccdf(s)
xl, yl = pow_linemaker(-0.5, (1, 0.6), 0.1, 10)

psl = scaling_plot(x, y; comparable=(xl, yl, "s ^ (1/2)"))
gui()
readline()