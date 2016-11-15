workspace()
using BoxChopper

box = sort(rand(3,2), 2);
r₀ = mean(box, 2)[:,1];
nout = randn(3);

@code_warntype BoxChopper.edgecepts(box, nout, r₀)
