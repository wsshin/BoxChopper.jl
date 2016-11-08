using BoxChopper
using Base.Test

const Float = typeof(0.0)

@testset "BoxChopper" begin
    @testset "edgecepts" begin
        box = [0 1; 0 1; 0 1]
        nout = [1,1,1]
        r₀ = [0,0,1]

        cepts, hascept_in, hascept_ex = BoxChopper.edgecepts(box, nout, r₀)

        cepts_expected = Array{Float,3}(3,2,2)
        for w = 1:3
            cepts_expected[w,:,:] = [1. 0.; 0. -1.]
        end
        @test cepts == cepts_expected

        hascept_in_expected = Array{Bool}(2,2,2,3)
        hascept_in_expected[1,1,1,:] = [true, true, true]
        hascept_in_expected[2,1,1,:] = [false, false, false]
        hascept_in_expected[1,2,1,:] = [false, false, false]
        hascept_in_expected[2,2,1,:] = [true, true, false]
        hascept_in_expected[1,1,2,:] = [false, false, false]
        hascept_in_expected[2,1,2,:] = [true, false, true]
        hascept_in_expected[1,2,2,:] = [false, true, true]
        hascept_in_expected[2,2,2,:] = [false, false, false]
        @test hascept_in == hascept_in_expected

        hascept_ex_expected = Array{Bool}(2,2,2,3)
        hascept_ex_expected = copy(hascept_in_expected)
        hascept_ex_expected[2,2,2,:] = [true, true, true]
        @test hascept_ex == hascept_ex_expected

        @test begin
            result = true
            for i = 1:100
                nout = rand(3)
                r₀ = (rand(3) + 1) / 2  # rand() ϵ [0, 1), so (rand()+1)/2 ϵ (0, 1), strictly within box
                ~, hascept_in, hascept_ex = BoxChopper.edgecepts(box, nout, r₀)

                Ncept_ex = squeeze(sum(hascept_ex, 4), 4)

                # Test the two vertices whose three edges cross the plane when extended are body-diagonally opposite.
                result &= all(map(sum, [findn(Ncept_ex .== 3)...]) .== 3)

                Ncept_in = squeeze(sum(hascept_in, 4), 4)
                t1, t2 = find(Ncept_ex .== 3)
                if Ncept_in[t1] < Ncept_in[t2]
                    (t1, t2) = (t2, t1)
                end
                N1, N2 = Ncept_in[[t1,t2]]  # edges from those body-diagonal vertices cross plane N1 and N2 times (without edge extension)
                Ntot = sum(hascept_in) / 2  # cross section is Ntot-gon.

                # Test the postulates about the relationship between the cross-sectional shape and (N1, N2).
                result &= (N2==0 && N1==0 && Ntot==6) || (N2==0 && N1==1 && Ntot==5) || (N2==0 && N1==2 && Ntot==4) || (N2==0 && N1==3 && Ntot==3) || (N2==1 && N1==1 && Ntot==4)
            end
            result
        end
    end  # @testset "edgecepts"

    @testset "contains" begin
        box = [0. 1.; 0. 1.; 0. 1.]
        nout = [1.,1.,1.]
        r₀ = [0.,0.,1.]
        @test BoxChopper.contains(nout, r₀, BoxChopper.vertices(box), true) == [true,true,true,false,true,false,false,false]
    end

    @testset "rvol_tricyl" begin
        box = [0 1; 0 1; 0 1]
        nout = [1, 1, 0]

        @test (r₀ = [0.5, 0., 0.]; volfrac(box, nout, r₀) ≈ 0.125)
        @test (r₀ = [1., 0., 0.]; volfrac(box, nout, r₀) ≈ 0.5)
        @test (r₀ = [1., 0., 0.]; nout = [1., 2., 0.]; volfrac(box, nout, r₀) ≈ 0.25)
    end

    @testset "rvol_quadcyl" begin
        @test begin
            result = true
            for i = 1:100
                box = sort(rand(3,2), 2)
                r₀ = mean(box, 2)[:,1]
                nout = randn(3)
                nout[rand(1:3)] = 0
                result &= volfrac(box, nout, r₀) ≈ 0.5
            end
            result
        end
    end

    @testset "volfrac" begin
        # Test random cases.
        @test begin
            result = true
            for i = 1:100
                box = sort(rand(3,2), 2)
                r₀ = mean(box, 2)[:,1]
                nout = randn(3)
                result &= volfrac(box, nout, r₀) ≈ 0.5
            end
            result
        end

        # Test boundary cases.
        box = [0 1; 0 1; 0 1]
        r₀ = [0,0,0]
        @test (nout = [1,0,0]; volfrac(box, nout, r₀) ≈ 0.)  # completely outside
        @test (nout = [-1,0,0]; volfrac(box, nout, r₀) ≈ 1.)  # completely inside
        @test (nout = [-1,1,0]; volfrac(box, nout, r₀) ≈ 0.5)  # rvol_tricyl()
        @test (nout = [-1,-1,1]; volfrac(box, nout, r₀) ≈ 5/6)  # rvol_gensect()
        @test (nout = [1,-2,1]; volfrac(box, nout, r₀) ≈ 0.5)  # rvol_quadsect()
        r₀ = [0.5,0.5,0]
        @test (nout = [-2,1,0]; volfrac(box, nout, r₀) ≈ 0.5)  # rvol_quadcyl()

        # Test rvol_quadsect() for nontrivial cases.
        box = [0 1; 0 1; 0 2]
        r₀ = [0.5, 0.5, 0.5]
        @test begin
            result = true
            for i = 1:100
                nout = [randn()/20, randn()/20, 1]
                result &= volfrac(box, nout, r₀) ≈ 0.5/2
            end
            result
        end
    end

end  # @testset "BoxChopper"
