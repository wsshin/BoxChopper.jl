# BoxChopper

When a plane runs across a 3-dimensional box, the box is divided into two pieces by the plane.  This Julia package calculates how large these two pieces are in volume with respect to the original box volume.

## Usage
The only function `BoxChopper` exports is `volfrac()`:
```
rvol = volfrac(box, nout, r₀, [verbose = false])
```

`box = [xn xp; yn yp; zn zp]` is a 3-by-2 matrix that indicates the locations of the x-, y-, z-boundaries of the box.  The box needs to be aligned with the Cartesian axes.

`nout = [nx, ny, nz]` and `r₀ = [rx, ry, rz]` are the normal vector of the plane and a point on the plane, respectively.  `nout` does not have to be normalized, but its sign has significance: among the two half-spaces divided by the plane, the one in the `nout` side of the plane is considered "outside", and the other is considered "inside".  

Now, among the two subpieces of the box divided by the plane, take the one in the "inside" portion of the space.  The return value `rvol` is the fraction of this "inside" piece's volume with respect to the original box's volume.  The volume fraction of the other piece is `1 - rvol`.  Both `rvol` and `1 - rvol` are between `0.0`  and `1.0`, inclusive.
