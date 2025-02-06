include("../src/JCommonCode.jl")

using .JCommonCode
using Plots

stress = readline_data("test/test_data/test_data_1.txt")
plot(stress; legend=false)
gui()
readline()

dstress = diff(stress)
pds = plot(dstress; legend=false)
gui()
readline()

v, t, s, d, b, e = get_slips(; vel=dstress, drops=true, threshold=-1)

for (av_vel, av_t) in zip(v, t) 
    
    av_vel .*= -1
    plot!(pds, av_t, av_vel; color=:red)
end
gui()
readline()
savefig(pds, "get_slips_vel.png")

x, y = ccdf(s)
xl, yl = pow_linemaker(-0.5, (1, 0.6), 0.1, 10)

psl = scaling_plot(x, y; comparable=(xl, yl, "s ^ (1/2)"))
gui()
readline()
savefig(psl, "CCDF_Plot_s0.5.png")