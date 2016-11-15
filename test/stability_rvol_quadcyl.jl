workspace()
using BoxChopper

box = [0 1; 0 1; 0 1]
r₀ = [0.5, 0.5, 0.5]
nout = [-2, 1, 0]

v = BoxChopper.vertices(box)
Nci = length(find(BoxChopper.contains(nout, r₀, v, true)))  # number of vertices contained in half-space, including boundary plane
Nnci = length(find(BoxChopper.contains(-nout, r₀, v, true)))  # number of vertices contained in half-space (≈ not contained (nc) in half space), including boundary plane
svc = BoxChopper.VI2VS[find(BoxChopper.contains(nout, r₀, v, false))]  # array of signs (sx,sy,sz) of vertices contained in half-space, excluding boundary plane
svnc = BoxChopper.VI2VS[find(BoxChopper.contains(-nout, r₀, v, false))]  # array of signs (sx,sy,sz) of vertices contained in the other half-space (≈ not contained (nc) in half space) excluding boundary plane
Nc = length(svc)  # number of vertices contained in half-space, excluding boundary plane
Nnc = length(svnc)  # number of vertices contained in the other half-space (≈ not contained (nc) in half space), excluding boundary plane

cepts, hascept_in, hascept_ex = BoxChopper.edgecepts(box, nout, r₀)

r = find(nout .== 0)[1]  # cylinder axis; if two components of nout are zero, choose first component's direction as cylinder axis

@code_warntype BoxChopper.rvol_quadcyl(box, cepts, hascept_in, r, svc)
