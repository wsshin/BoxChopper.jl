workspace();
using BoxChopper;

box = sort(randn(3,2), 2)
r₀ = mean(box, 2)[:,1]
nout = randn(3)

@time volfrac(box, nout, r₀)
