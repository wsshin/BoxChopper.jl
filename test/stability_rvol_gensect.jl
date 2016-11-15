workspace();
using BoxChopper;

box = [0 1; 0 1; 0 1]
r₀ = [0,0,0]
nout = [-1,-1,1]

v = BoxChopper.vertices(box)
Nci = length(find(BoxChopper.contains(nout, r₀, v, true)))  # number of vertices contained in half-space, including boundary plane
Nnci = length(find(BoxChopper.contains(-nout, r₀, v, true)))  # number of vertices contained in half-space (≈ not contained (nc) in half space), including boundary plane
svc = BoxChopper.VI2VS[find(BoxChopper.contains(nout, r₀, v, false))]  # array of signs (sx,sy,sz) of vertices contained in half-space, excluding boundary plane
svnc = BoxChopper.VI2VS[find(BoxChopper.contains(-nout, r₀, v, false))]  # array of signs (sx,sy,sz) of vertices contained in the other half-space (≈ not contained (nc) in half space) excluding boundary plane
Nc = length(svc)  # number of vertices contained in half-space, excluding boundary plane
Nnc = length(svnc)  # number of vertices contained in the other half-space (≈ not contained (nc) in half space), excluding boundary plane

cepts, hascept_in, hascept_ex = BoxChopper.edgecepts(box, nout, r₀)

Ncept_ex = sum(hascept_ex, 4)[:,:,:,1]  # number of each vertex's edges that cross plane when extended
indc = find(Ncept_ex .== 3)
sv1, sv2 = BoxChopper.VI2VS[indc[1]], BoxChopper.VI2VS[indc[2]]

Ncept_in = sum(hascept_in, 4)[:,:,:,1]  # number of each vertex's edges that cross plane
N1 = Ncept_in[sv1...]
N2 = Ncept_in[sv2...]
if N1 < N2
    sv1, sv2 = sv2, sv1
    N1, N2 = N2, N1
end

v1 = [box[w,sv1[w]] for w = BoxChopper.AXES]
@code_warntype BoxChopper.rvol_gensect(box, cepts, hascept_in, sv1)
