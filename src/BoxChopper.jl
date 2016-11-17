module BoxChopper

export volfrac  # function

const Float = typeof(0.)

const N, P = 1, 2  # negative, positive directions
const NP = [N,P]

const X, Y, Z = 1, 2, 3
const XYZ, YZX, ZXY = [X,Y,Z], [Y,Z,X], [Z,X,Y]  # array is faster to construct than tuple
const CYC_AXES = [XYZ, YZX, ZXY]

typealias AbsVec AbstractVector
typealias AbsMat AbstractMatrix

function calc_vcbits{T<:Real,S<:Real}(box::AbsMat{T}, nout::AbsVec{S}, nr₀::Real)
    vcbits = 0x00
    bit = 0x01
    for sz = NP
        for sy = NP
            for sx = NP
                if nout[X]*box[X,sx] + nout[Y]*box[Y,sy] + nout[Z]*box[Z,sz] ≤ nr₀
                    vcbits |= bit
                end
                bit <<= 1
            end
        end
    end

    return vcbits
end

function isquadsect(vcbits::UInt8)
    return vcbits==0x0F || vcbits==~0x0F || vcbits==0x33 || vcbits==~0x33 || vcbits==0x55 || vcbits==~0x55
    # n = count_ones(vcbits)
    # m = count_ones(vcbits & 0x0F)
    # return n==4 && iseven(m)
end

function edgedir_quadsect(vcbits::UInt8)
    if vcbits==0x0F || vcbits==~0x0F
        return Z
    elseif vcbits==0x33 || vcbits==~0x33
        return Y
    elseif vcbits==0x55 || vcbits==~0x55
        return X
    else
        throw(ArgumentError("vcbits = $(bits(vcbits)) should correspond to quadrangular cross section."))
    end
end

function rvol_quadsect{T<:Real,S<:Real}(box::AbsMat{T}, nout::AbsVec{S}, nr₀::Real, vcbits::UInt8)
    # Return the volume fraction for the case where the plane croses four parallel edges.
    # w: direction of four parallel edges that cross plane

    w = edgedir_quadsect(vcbits)
    ~, u, v = CYC_AXES[w]
    ∆w = box[w,P] - box[w,N]
    mean_cepts = 4nr₀
    for sv = NP
        for su = NP
            mean_cepts -= nout[u]*box[u,su] + nout[v]*box[v,sv]
        end
    end
    mean_cepts /=  nout[w]*4∆w

    sw = nout[w]>0 ? N:P

    return abs(mean_cepts - box[w,sw]/∆w)
end

function rvol_gensect{T<:Real,S<:Real}(box::AbsMat{T}, nout::AbsVec{S}, nr₀::Real, vcbits::UInt8)
    # Calculate the volume fraction of most general cases, by cutting out corners
    # from a triangular pyramid.
    # Assume count_ones(vcbits) ≤ 4.  Othewise, call this function with flipped
    # nout, nr₀, vcbits.

    s = [nout[X]≥0 ? N:P, nout[Y]≥0 ? N:P, nout[Z]≥0 ? N:P]  # signs of corner
    c = [box[X,s[X]], box[Y,s[Y]], box[Z,s[Z]]]  # corner coordinates
    ∆ = [box[X,P]-box[X,N], box[Y,P]-box[Y,N], box[Z,P]-box[Z,N]]  # box edges
    nc = [nout[X]*c[X], nout[Y]*c[Y], nout[Z]*c[Z]]  # nout * corner
    r = [abs((nr₀-nc[Y]-nc[Z])/nout[X]-c[X])/∆[X], abs((nr₀-nc[Z]-nc[X])/nout[Y]-c[Y])/∆[Y], abs((nr₀-nc[X]-nc[Y])/nout[Z]-c[Z])/∆[Z]]

    sort!(r)
    r3_ = 1 - 1/r[3]
    rvol_core = 1+r3_+r3_^2
    if r[2] > 1
        r2_ = 1 - 1/r[2]
        rvol_core -= r[3]*r2_^3
    end
    if r[1] > 1
        r1_ = 1 - 1/r[1]
        rvol_core -= r[3]*r1_^3
    end

    return rvol_core * r[1] * r[2] / 6
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

    # If box is completely within the half-space, rvol = 1.
    # if box is completely outside the half-space, rvol = 0.
    if nvc == 8; return 1.; end
    if nvc == 0; return 0.; end

    if isquadsect(vcbits); return rvol_quadsect(box, nout, nr₀, vcbits); end

    needflip = false
    if nvc ≥ 5
        needflip = true
        nout = -nout
        nr₀ = -nr₀
        vcbits = ~vcbits
        nvc = 8-nvc
    end

    return needflip ? 1 - rvol_gensect(box, nout, nr₀, vcbits) : rvol_gensect(box, nout, nr₀, vcbits)
end

end  # module PlanarBoxDivider
