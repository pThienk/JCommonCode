include("../src/JCommonCode.jl")

using .JCommonCode
using Plots

stress = readline_data("test/test_data/test_data_2.txt")
dstress = diff(stress)

v, t, s, d, b, e = get_slips(; vel=dstress, drops=true, threshold=-1)

x, y = ccdf(s)
#xl, yl = pow_linemaker(-0.5, (1, 0.6), 0.1, 10)

#psl = scaling_plot(x, y; comparable=(xl, yl, "s ^ (1/2)"))

allt, allv, avgt, avgv, avgstd, centers, width = shapes(v, t, s, d; style=:duration_shapes, nbins=5)

println("width=$width, " * "centers=$(join(string.(centers), ","))")

plt = plot(avgt[3], avgv[3]; title="Scaling Collapse - Normal Weakening - e=0.01", color=:auto, ls=:auto, lw=2, label="T-center=$(Float16(centers[3]))",
 legend=:outertop, legendcolumns=2)

plot!(avgt[4], avgv[4]; label="T-center=$(Float16(centers[4]))", color=:auto, ls=:auto, lw=2)
plot!(avgt[5], avgv[5]; label="T-center=$(Float16(centers[5]))", color=:auto, ls=:auto, lw=2)

gui()
readline()