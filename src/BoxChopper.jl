module BoxChopper

export volfrac  # function

const Float = typeof(0.0)
const Nv = 8  # number of vertices of box
const Ne = 12  # number of edges of box
const X, Y, Z = [1, 2, 3]
const N, P = [1, 2]  # negative, positive directions
const AXES = [X, Y, Z]
const SIGNS = [N, P]

# If VE2E[sx,sy,sz,w] = [w,su,sv], then the edge corresponding to hascept[sx,sy,sz,w]
# is the same as the edge corresponding to cepts[w,su,sv].
const VE2E = [(
    s = [sx,sy,sz];
    [w, s[w%3+1], s[(w+1)%3+1]]
) for sx = SIGNS, sy = SIGNS, sz = SIGNS, w = AXES]

const VI2VS = [[ind2sub((2,2,2),li)...] for li = 1:8]  # vertex index to subscripts

function floatize_box{T<:Real}(box::Array{T,2})
    # box: [xn xp; yn yp; zn zp] (boundary locations of box)
    size(box) ≠ (3,2) && throw(ArgumentError("box = $box should be size-(3,3)."))
    return T==Float ? box : float(box)
end

function floatize_plane{T<:Real,S<:Real}(nout::Array{T,1}, r₀::Array{S,1})
    # nout: [nx, ny, nz] (outward normal of plane)
    # r₀: [x₀, y₀, z₀] (point on plane)
    # plane equation: p(r) = nout' * (r - r₀) = 0 (p > 0: nout portion, p < 0: -nout portion)
    length(nout) ≠ 3 && throw(ArgumentError("nout = $nout should be length-3."))
    all(nout .== 0) && throw(ArgumentError("nout = $nout should be nonzero vector."))
    length(r₀) ≠ 3 && throw(ArgumentError("r₀ = $r₀ should be length-3."))

    return (T==Float ? nout:float(nout), S==Float ? r₀:float(r₀))
end

function edgecepts(box, nout, r₀)
    box = floatize_box(box)
    nout, r₀ = floatize_plane(nout, r₀)

    # 8 vertices are listed in the order where x is the fastest varying coordinate,
    # y is the next, and z is the slowest.  In other words,
    # [xn, yn, zn], [xp, yn, zn], [xn, yp, zn], [xp, yp, zn],
    # [xn, yn, zp], [xp, yn, zp], [xn, yp, zp], [xp, yp, zp].

    # 12 edges are listed in the order of four x-directional edges, then four y-
    # directional edges, and then four z-directional edges.  In each group, the
    # four edges are stored in a 2×2 array, whose faster varying index is the coordinate
    # that cyclically comes next to the edge direction.  (E.g., for the y-directional
    # edges, the first index varies the z-coordinate, and the second index varies
    # the x-coordinate.)

    # The function first returns 12 intercepts between 12 edges and the plane.
    # Each intercept is described by a single coordinate, because the box is aligned
    # with the Cartesian axes and therefore the other two coordinates are trivial
    # from the box.

    # The function returns two additional 2-by-2-by-2-by-3 arrays of Boolean values.
    # The first three indices specify a vertex, and the last index specifies edge
    # direction.  (E.g., indices [N,N,N,X] specifies the x-directional edge originating
    # from the vertex at [xn,yn,zn].)

    # In the first array, the Boolean entries indicate whether the edge specified
    # by the indices crosses the plane.  In the second array, the Boolean entries
    # indicate whether the "ray" (or a half line) emanating from the vertex in the
    # edge direction crosses the plane.

    # For both arrays, the edges or rays exclude the vertices from which they originate.
    # For example, for the edge from [0,0,0] to [1,0,0], the points [0.5,0,0] and
    # [1,0,0] are contained in the edge, but [0,0,0] is considered not.

    nr₀ = nout⋅r₀

    # Below, note that the plane equation is nout[X]*x + nout[Y]*y + nout[Z]*z - nout⋅r₀ = 0
    cepts = [(
        (u, v) = (w%3+1, (w+1)%3+1);
        ruv = [box[u,su], box[v,sv]];
        nout[w]==0.0 ? NaN:(nr₀ - nout[[u,v]]⋅ruv) / nout[w]
    ) for w = AXES, su = SIGNS, sv = SIGNS]  # about assignment of NaN, see below

    # hascept_ex = Array{Bool,4}(2,2,2,3)  # does edge cross plane when extended?
    # hascept_in = Array{Bool,4}(2,2,2,3)  # does edge cross plane?
    hascept_ex = Array{Bool}(2,2,2,3)  # does edge cross plane when extended?
    hascept_in = Array{Bool}(2,2,2,3)  # does edge cross plane?

    # Comparison (<, >, ≤, ≥) of any number with NaN is always false.  Therefore,
    # if an edge of the box is parallel to the plane, its hascept is false (i.e.,
    # does not cross the plane).
    exop = (<, >)  # exclusive inequality operators
    inop = (≤, ≥)  # inclusive inequality operators

    hascept_ex = [(
        s = [sx,sy,sz];
        (u, v) = (w%3+1, (w+1)%3+1);
        exop[s[w]](box[w,s[w]], cepts[w,s[u],s[v]])
    ) for sx = SIGNS, sy = SIGNS, sz = SIGNS, w = AXES]  # does edge cross plane when extended?

    hascept_in = [(
        s = [sx,sy,sz];
        (u, v) = (w%3+1, (w+1)%3+1);
        negsw = s[w]%2+1;
        hascept_ex[sx,sy,sz,w] && inop[s[w]](cepts[w,s[u],s[v]], box[w,negsw])
    ) for sx = SIGNS, sy = SIGNS, sz = SIGNS, w = AXES]  # does edge cross plane?

    # Need to return a Boolean array of size = size(cepts).
    # Or, some data structure that connects hascept to the indices of cepts (i.e., edge indices)?
    # When a vertex index and its x, y, or z is given, we should be able to get
    # the edge index used to index cepts.  This can be a global parameter, because
    # it is problem-independent.  The mapping should be of size = size(hascept_in)
    # and Int type.
    return cepts, hascept_in, hascept_ex
end

function vertices(box)
    # Return the coordinates of the vertices of the given box.

    box = floatize_box(box)
    return [box[X,N] box[Y,N] box[Z,N];
            box[X,P] box[Y,N] box[Z,N];
            box[X,N] box[Y,P] box[Z,N];
            box[X,P] box[Y,P] box[Z,N];
            box[X,N] box[Y,N] box[Z,P];
            box[X,P] box[Y,N] box[Z,P];
            box[X,N] box[Y,P] box[Z,P];
            box[X,P] box[Y,P] box[Z,P]]
end

function contains(nout::Array{Float,1}, r₀::Array{Float,1}, pt::Array{Float,2}, isinclusive)
    # Tells if pt is contained in the half-space whose outward normal is nout.
    size(pt)[2] ≠ 3 && throw(ArgumentError("pt = $pt should have three columns."))
    return (isinclusive ? (.≤) : (.<))(broadcast(-, pt, transpose(r₀)) * nout, 0.)
end

function contains(nout::Array{Float,1}, r₀::Array{Float,1}, pt::Array{Float,1}, isinclusive)
    length(pt) ≠ 3 && throw(ArgumentError("pt = $pt should be length-3."))
    return (isinclusive ? (≤) : (<))((pt-r₀)⋅nout, 0.)
end

function rvol_tricyl(box, cepts, hascept_in, r::Int, sv::Array{Array{Int,1},1})
    # Return the volume fraction of the triangular cylinder in the box.
    # r: cylinder axis (one of X, Y, Z)
    # sv: array of signs [sx,sy,sz] of vertices

    length(sv) ≠ 2 && throw(ArgumentError("sv = $sv should be length-2."))

    sv1, sv2 = sv

    axes = [r%3+1, (r+1)%3+1]  # axes normal to cylinder axis
    !all(hascept_in[sv1...,axes]) && throw(ArgumentError("All edges in $axes-direction from vertex indexed by $sv1 should cross plane."))
    !all(hascept_in[sv2...,axes]) && throw(ArgumentError("All edges in $axes-direction from vertex indexed by $sv2 should cross plane."))

    eind1 = VE2E[sv1...,axes]  # indices of r-normal edges from vertex sv1
    eind2 = VE2E[sv2...,axes]  # indices of r-normal edges from vertex sv2
    cepts[eind1[1]...] ≠ cepts[eind2[1]...] && throw(ArgumentError("Edges in $(axes[1])-direction from vertices indexed by $sv1 and $sv2 should have same intercept location."))
    cepts[eind1[2]...] ≠ cepts[eind2[2]...] && throw(ArgumentError("Edges in $(axes[2])-direction from vertices indexed by $sv1 and $sv2 should have same intercept location."))

    ∆w = (box[:,2]-box[:,1])[axes]
    d = [abs(cepts[eind1[w]...] - box[w,sv1[w]]) for w = axes]

    return prod(d./∆w)/2.
end

function rvol_quadcyl(box, cepts, hascept_in, r::Int, sv::Array{Array{Int,1},1})
    # Return the volume fraction of the quadrangular cylinder in the box.
    # r: cylinder axis (one of X, Y, Z)
    # sv: array of signs [sx,sy,sz] of vertices

    length(sv) ≠ 4 && throw(ArgumentError("sv = $sv should be length-4."))

    p = find(sv[1].==sv[2].==sv[3].==sv[4])  # direction normal to plane of four vertices
    length(p) ≠ 1 && throw(ArgumentError("Vertices indexed by sv = $sv should be on same face."))

    p = p[1]
    p == r && throw(ArgumentError("Face containing vertices indexed by sv = $sv should be parallel to r = $r-direction."))

    q = p%3+1
    if q==r; q = q%3+1; end  # direction (other than r) in plane of four vertices

    sv1 = sv[1]  # pick any one vertex contained in half space
    sv2 = copy(sv[1])
    sv2[q] = 3 - sv1[q]  # vertex in q-directino from sv1  (if sv1[q] == N, then sv2[q] == P, and vice versa)
    length(find(sv.==[sv2])) ≠ 1 && throw(ArgumentError("sv = $sv should contain $sv2 once and only once."))

    ep1 = VE2E[sv1..., p]
    ep2 = VE2E[sv2..., p]

    ∆p = box[p,2] - box[p,1]
    a1 = abs(cepts[ep1...] - box[p,sv1[p]])
    a2 = abs(cepts[ep2...] - box[p,sv2[p]])

    return (a1+a2)/2∆p
end

sg2(x) = 1 .+ x .+ x.^2

function get_rs(box, cepts, sv1::Array{Int,1})
    eind = VE2E[sv1...,:]  # indices of three edges from vertex sv1

    ∆w = box[:,2] - box[:,1]
    d = [abs(cepts[eind[w]...] - box[w,sv1[w]]) for w = AXES]
    r = d ./ ∆w
    r_ = 1. - ∆w./d
    sg2r_ = sg2(r_)

    return (r, sg2r_)
end

function rvol_gensect(box, cepts, hascept_in, sv1::Array{Int,1})
    # Calculate the volume fraction of most general cases, by cutting out corners
    # from a triangular pyramid.
    # sv1: signs [sx,sy,sz] of vertex

    r, sg2r_ = get_rs(box, cepts, sv1)

    icr = hascept_in[sv1...,:]  # is this edge crossing plane?
    xncr = AXES[!icr]  # axes not crossing plane
    Nncr = length(xncr)  # number of edges not crossing plane
    rvol = 0.
    for w = xncr
        u, v = w%3+1, (w+1)%3+1  # complementary of w
        rvol += r[u]*r[v]*sg2r_[w] / 6.
    end

    rvol += (1-Nncr) * prod(r) / 6.

    return rvol
end

function rvol_quadsect(box, cepts, hascept_in, sv1::Array{Int,1})
    # Return the volume fraction for the case where the plane croses four parallel edges.
    # sv1: signs [sx,sy,sz] of vertex

    icr = hascept_in[sv1...,:]  # is this edge crossing plane?
    sum(icr) ≠ 1 && throw(ArgumentError("One and only one edge from vertex indexed by sv1 = $sv1 should cross plane."))

    r = AXES[icr][1]  # axis normal to face that do not cross plane
    ∆r = box[r,2] - box[r,1]

    return sum(abs(cepts[r,:,:] - box[r,sv1[r]])) / 4∆r
end

"""
    volfrac(box, nout, r₀, [verbose = false])

Returns the volume fraction `rvol = vol(box ⋂ half-space) / vol(box)` that indicates
how much portion of the volume of a box is included the half-space defined by a
plane.  The result is between 0.0 and 1.0, inclusive.

The box needs to be aligned with the Cartesian axes, and is described by `box =
[xn xp; yn yp; zn zp]` that specifies the x-, y-, z-boundaries of the box.  The
half-space is described by the boundary plane.  The boundary plane is described by
the outward normal vector `nout = [nx, ny, nz]` and a point `r₀ = [rx, ry, rz]`
on the plane.  `nout` does not have to be normalized.

`volfrac()` uses different volume calculation functions internally depending on
the relative positions and directions of the box and the plane.  The optional ArgumentError
`verbose` prints out which of these functions was used.  The default value of `verbose`
is `false`, meaning that it does not print out the information.
"""
function volfrac(box, nout, r₀, verbose = false)
    # Return the volume fraction rvol = vol(box ⋂ half-space) / vol(box).

    box = floatize_box(box)
    nout, r₀ = floatize_plane(nout, r₀)

    # If box is completely within the half-space, rvol = 1.
    # if box is completely outside the half-space, rvol = 0.
    v = vertices(box)
    Nci = length(find(contains(nout, r₀, v, true)))  # number of vertices contained in half-space, including boundary plane
    Nnci = length(find(contains(-nout, r₀, v, true)))  # number of vertices contained in half-space (≈ not contained (nc) in half space), including boundary plane
    svc = VI2VS[find(contains(nout, r₀, v, false))]  # array of signs [sx,sy,sz] of vertices contained in half-space, excluding boundary plane
    svnc = VI2VS[find(contains(-nout, r₀, v, false))]  # array of signs [sx,sy,sz] of vertices contained in the other half-space (≈ not contained (nc) in half space) excluding boundary plane
    Nc = length(svc)  # number of vertices contained in half-space, excluding boundary plane
    Nnc = length(svnc)  # number of vertices contained in the other half-space (≈ not contained (nc) in half space), excluding boundary plane

    # Note that Nnci is not necessarily 8 - Nci, because the two half-spaces considered
    # to cound Nci and Nnci are not mutually exclusive: they share the plane.
    if Nci == 8  # if all vertices are contained in half-space, then rvol = 1
        if verbose println("Completely inside.") end
        return 1.
    end
    if Nnci == 8  # if all vertices are contained in the other half-space, then rvol = 0
        if verbose println("Completely outside.") end
        return 0.
    end  # after this, 0 < rvol < 1.

    cepts, hascept_in, hascept_ex = edgecepts(box, nout, r₀)

    # If nout is in the xy-, yz-, or zx-plane, box ⋂ half-space is a cylinder.
    if any(nout .== 0)
        r = find(nout .== 0)[1]  # cylinder axis; if two components of nout are zero, choose first component's direction as cylinder axis

        if Nc == 2 || Nnc == 2  # cylinder base is triangle
            if Nc == 2
                sv = svc
                needflip = false
            else  # Nnc == 2
                sv = svnc
                needflip = true
            end
            if verbose println("Use rvol_tricyl.") end
            rvol = rvol_tricyl(box, cepts, hascept_in, r, sv)
            return needflip ? 1.-rvol : rvol
        else  # cylinder base is trapezoid
            assert(Nc == Nnc == 4)
            if verbose println("Use rvol_quadcyl.") end
            return rvol_quadcyl(box, cepts, hascept_in, r, svc)
        end
    end

    # Below, deal with all the other general cases.

    # Find two vertices whose all three edges, when extended in the direction
    # from the vertex, cross the plane.
    Ncept_ex = sum(hascept_ex, 4)[:,:,:,1]  # number of each vertex's edges that cross plane when extended
    sv1, sv2 = VI2VS[find(Ncept_ex .== 3)]  # signs [sx,sy,sz] of vertices whose edges cross plane three times when extended
    assert(all(sv1+sv2 .== 3))  # sv1 and sv2 are body-diagonally opposite

    # Between sv1 and sv2, make sv1 the vertex with more edges crossing the
    # plane without extension.
    Ncept_in = sum(hascept_in, 4)[:,:,:,1]  # number of each vertex's edges that cross plane
    N1 = Ncept_in[sv1...]
    N2 = Ncept_in[sv2...]
    if N1 < N2
        sv1, sv2 = sv2, sv1
        N1, N2 = N2, N1
    end

    # Determine the shape of the plane-box intersection from N1 and N2, and calculate
    # rvol for each cross-sectional shape accordingly.

    # It is tempting to determine the shape of the cross-sectional polygon by simply
    # counting the number of edges crossing the plane, but there are a couple of
    # problems in this method.
    # 1) There are two kinds of quadrangular cross sections (one crossing four parallel
    # edges and the other crossing two pairs of two joinging edges), and simply
    # counting the number of plane-crossing edges cannot distinguish the two cases.
    # Of course, we can distinguish the two by checking whether the plane-crossing
    # edges are all parallel or not, but there is a way that does not require this
    # extra checking step.
    # 2) Determination of the cross-sectional shape is not the end of the story,
    # and we need to proceed to calculate the valume.  So, we should think about
    # what information facillitates volume calculation.  Most different kinds of
    # volume turn out to reduce to a right-angled triangular pyramid whose corners
    # are cut off.  Therefore, it is crucial to find the vertex that forms the right-angled
    # corner of the pyramid.  It turns out that this vertex information can also
    # be used to distinguish the cross-sectional shapes, including the two kinds
    # of the quadrangular cross sections.  If the vertex where the right-angled
    # corner of the pyramid is formed is needed later anyway and it can be also
    # used to determine the cross-sectional shape, there is no reason not to determine
    # the cross sectional shape by using this vertex information instead of counting
    # the number of edges crossing the plane!

    v1 = [box[w,sv1[w]] for w = AXES]
    assert((v1-r₀)⋅nout ≠ 0.)
    needflip = !contains(nout, r₀, v1, false)
    if N2 == 0  # cross section is hexagon | pentagon | quadrangle (with plane crossing two pairs of joining edges) || triangle
        assert(N1==0 || N1==1 || N1==2 || N1==3)
        if verbose println("Use rvol_gensect.") end
        rvol = rvol_gensect(box, cepts, hascept_in, sv1)
    else  # cross section is quadrangle (with plane crossing four parallel edges)
        assert(N2==1 && N1==1)
        if verbose println("Use rvol_quadsect.") end
        rvol = rvol_quadsect(box, cepts, hascept_in, sv1)
    end

    return needflip ? 1.-rvol : rvol
end

end  # module PlanarBoxDivider
