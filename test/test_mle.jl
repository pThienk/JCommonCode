include("../src/JCommonCode.jl")

using .JCommonCode
using Random

seed = 314159
x0 = 1
α = 0.5

Random.seed!(seed)
sample::Vector = [x0*rand()^(-1/α) for _ in 1:10e6]

stat = mle(sample)
print("ahat: " * string(stat[1]) * ", stderr: " * string(stat[2]) * "\n")

powers::Vector = scan_ks(sample, 0.4:0.01:0.6)
print("likely powers: " * join(powers, ", ") * "\n")

result = scan_mle(sample)
show(result)