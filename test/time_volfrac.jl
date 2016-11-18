workspace();
using BoxChopper;

nbox = 100^2

boxes = [sort(randn(3,2),2) for i = 1:nbox]
r₀s = map(box->mean(box,2)[:,1], boxes)
nouts = [randn(3) for i = 1:nbox]

volfrac(boxes[1], nouts[1], r₀s[1])  # execute this to initiate compilation

@time begin
    for i = 1:nbox
        volfrac(boxes[i], nouts[i], r₀s[i])
    end
end
