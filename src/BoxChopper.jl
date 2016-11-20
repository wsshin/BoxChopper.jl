module BoxChopper

export volfrac  # function

const N, P = 1, 2  # negative, positive directions
const NP = (N,P)  # tuple allocation is more efficient than array allocation

const X, Y, Z = 1, 2, 3  # x, y, z coordinates
const XYZ, YZX, ZXY = (X,Y,Z), (Y,Z,X), (Z,X,Y)
const CYC_AXES = (XYZ, YZX, ZXY)

typealias AbsVec AbstractVector
typealias AbsMat AbstractMatrix

function calc_vcbits{T<:Real,S<:Real}(box::AbsMat{T}, nout::AbsVec{S}, nr₀::Real)
    # Calculate the bit array that indicates vertex contained-ness.

    # The bit array is 8 bits (= 1 byte).  The kth bit is 1 if the kth verex is
    # contained in the half space defined by the plane, an 0 otherwise.  The 8 vertices
    # of the box are indexed 1 through 8 in the order of (---), (+--), (-+-), (++-),
    # (--+), (+-+), (-++), (+++), where, e.g., (+--) indicates the vertex at the
    # +x, -y, -z corner.
    nx, ny, nz = nout
    vcbits = 0x00
    bit = 0x01
    for sz = NP, sy = NP, sx = NP
        if nx*box[X,sx] + ny*box[Y,sy] + nz*box[Z,sz] ≤ nr₀
            vcbits |= bit
        end
        bit <<= 1
    end

    return vcbits
end

function isquadsect(vcbits::UInt8)
    # Determine if the vertex contained-ness bit array corresponds to the case where
    # the plane crosses one of the three sets of four parallel edges of the box.

    return vcbits==0x0F || vcbits==~0x0F || vcbits==0x33 || vcbits==~0x33 || vcbits==0x55 || vcbits==~0x55

    # Equivalent to
    # n = count_ones(vcbits)
    # m = count_ones(vcbits & 0x0F)
    # return n==4 && iseven(m)
end

function edgedir_quadsect(vcbits::UInt8)
    # For the cases where the plane crosses a set of four parallel edges of the
    # box, determine which direction those edges lie.

    if vcbits==0x0F || vcbits==~0x0F
        dir = Z
    elseif vcbits==0x33 || vcbits==~0x33
        dir = Y
    else
        assert(vcbits==0x55 || vcbits==~0x55)
        dir = X
    end

    return dir
end

function rvol_quadsect{T<:Real,S<:Real}(box::AbsMat{T}, nout::AbsVec{S}, nr₀::Real, vcbits::UInt8)
    # Return the volume fraction for the case where the plane croses a set of four
    # parallel edges.

    w = edgedir_quadsect(vcbits)
    ∆w = box[w,P] - box[w,N]

    ~, u, v = CYC_AXES[w]
    nw, nu, nv = nout[w], nout[u], nout[v]
    mean_cepts = 4nr₀
    for sv = NP, su = NP
        mean_cepts -= nu*box[u,su] + nv*box[v,sv]
    end
    mean_cepts /=  nw * 4∆w

    sw = nw>0 ? N:P  # nw cannot be 0
    return abs(mean_cepts - box[w,sw]/∆w)
end

function rvol_gensect{T<:Real,S<:Real}(box::AbsMat{T}, nout::AbsVec{S}, nr₀::Real, vcbits::UInt8)
    # Calculate the volume fraction of most general cases, by cutting out corners
    # from a triangular pyramid.
    # Assume count_ones(vcbits) ≤ 4.  Othewise, call this function with flipped
    # nout, nr₀, vcbits.

    nx, ny, nz = nout
    sx, sy, sz = (nx≥0 ? N:P), (ny≥0 ? N:P), (nz≥0 ? N:P)  # signs of corner
    cx, cy, cz = box[X,sx], box[Y,sy], box[Z,sz]  # corner coordinates
    ∆x, ∆y, ∆z = box[X,P]-box[X,N], box[Y,P]-box[Y,N], box[Z,P]-box[Z,N]  # box edges
    nxcx, nycy, nzcz = nx*cx, ny*cy, nz*cz

    rmax, rmid, rmin =  # (lengths from corner to intercetps) / (box edges)
        abs((nr₀-nycy-nzcz)/nx-cx)/∆x, abs((nr₀-nzcz-nxcx)/ny-cy)/∆y, abs((nr₀-nxcx-nycy)/nz-cz)/∆z

    # Sort rmax, rmin, rmin properly.
    if rmax < rmid; rmax, rmid = rmid, rmax; end
    if rmid < rmin; rmid, rmin = rmin, rmid; end
    if rmax < rmid; rmax, rmid = rmid, rmax; end

    # Calculate the volume of the triangular pyramid, and cut off appropriate corners.
    tmax = 1 - 1/rmax
    rvol_core = 1 + tmax + tmax^2
    if rmid > 1
        tmid = 1 - 1/rmid
        rvol_core -= rmax * tmid^3
    end
    if rmin > 1
        tmin = 1 - 1/rmin
        rvol_core -= rmax * tmin^3
    end

    return rvol_core * rmin * rmid / 6
end

"""
    volfrac(box, nout, r₀)

Returns the volume fraction `rvol = vol(box ⋂ half-space) / vol(box)` that indicates
how much portion of the volume of a box is included the half-space defined by a
plane.  The result is between 0.0 and 1.0, inclusive.

The box needs to be aligned with the Cartesian axes, and is described by `box =
[xn xp; yn yp; zn zp]` that specifies the x-, y-, z-boundaries of the box.  The
half-space is described by the boundary plane.  The boundary plane is described by
the outward normal vector `nout = [nx, ny, nz]` and a point `r₀ = [rx, ry, rz]`
on the plane.  `nout` does not have to be normalized.
"""
function volfrac{T<:Real,S<:Real,R<:Real}(box::AbsMat{T}, nout::AbsVec{S}, r₀::AbsVec{R})
    # Return the volume fraction rvol = vol(box ⋂ half-space) / vol(box).

    nr₀ = nout⋅r₀
    vcbits = calc_vcbits(box, nout, nr₀)
    nvc = count_ones(vcbits)  # number of vertices contained

    if nvc == 8  # box is inside half-space
        rvol = 1.
    elseif nvc == 0  # box is outside half-space
        rvol = 0.
    elseif isquadsect(vcbits) # plane crosses a set of four parallel edges of box
        rvol = rvol_quadsect(box, nout, nr₀, vcbits)
    elseif nvc ≤ 4 # general cases with nvc ≤ 4
        rvol = rvol_gensect(box, nout, nr₀, vcbits)
    else  # general cases with nvc ≥ 5
        assert(nvc ≥ 5)
        rvol = 1. - rvol_gensect(box, -nout, -nr₀, ~vcbits)
    end

    return rvol
end

end  # module PlanarBoxDivider
