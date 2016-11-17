workspace()
using BoxChopper

box = [0 1; 0 1; 0 1]
r₀ = [0.5, 0.5, 0.5]
nout = [1, -2, 1]

@code_warntype volfrac(box, nout, r₀)
