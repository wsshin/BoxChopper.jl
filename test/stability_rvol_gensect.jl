workspace();
using BoxChopper;

box = [0 1; 0 1; 0 1]
r₀ = [0,0,0]
nout = [-1,-1,1]

nr₀ = nout⋅r₀
vcbits = BoxChopper.calc_vcbits(box, nout, nr₀)
nvc = count_ones(vcbits)  # number of vertices contained

needflip = false
if nvc ≥ 5
    needflip = true
    nout = -nout
    nr₀ = -nr₀
    vcbits = ~vcbits
    nvc = 8-nvc
end

@code_warntype BoxChopper.rvol_gensect(box, nout, nr₀, vcbits)
