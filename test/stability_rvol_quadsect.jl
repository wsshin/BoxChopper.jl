workspace();
using BoxChopper;

box = [0 1; 0 1; 0 1]
r₀ = [0.5, 0.5, 0.5]
nout = [-2,-1,0]

nr₀ = nout⋅r₀
vcbits = BoxChopper.calc_vcbits(box, nout, nr₀)
nvc = count_ones(vcbits)  # number of vertices contained

assert(BoxChopper.isquadsect(vcbits))
@code_warntype BoxChopper.rvol_quadsect(box, nout, nr₀, vcbits)
