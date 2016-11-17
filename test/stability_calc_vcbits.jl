workspace();
using BoxChopper;

box = [0 1; 0 1; 0 1]
r₀ = [0,0,0]
nout = [-1,-1,1]

nr₀ = nout⋅r₀
@code_warntype calc_vcbits(box, nout, nr₀)
