workspace();
using BoxChopper;

box = sort(randn(3,2), 2)
r₀ = mean(box, 2)[:,1]
nout = randn(3)

volfrac(box, nout, r₀)  # execute this to initiate compilation

@time begin
    for i = 1:100^3
        volfrac(box, nout, r₀)
    end
end
