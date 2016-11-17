using BoxChopper
using Base.Test

const Float = typeof(0.0)

@testset "BoxChopper" begin
    @testset "triangular cylinder" begin
        box = [0 1; 0 1; 0 1]
        nout = [1, 1, 0]

        @test (r₀ = [0.5, 0., 0.]; volfrac(box, nout, r₀) ≈ 0.125)
        @test (r₀ = [0.5, 0., 0.]; volfrac(box, -nout, r₀) ≈ 0.875)
        @test (r₀ = [1., 0., 0.]; volfrac(box, nout, r₀) ≈ 0.5)
        @test (r₀ = [1., 0., 0.]; volfrac(box, -nout, r₀) ≈ 0.5)
        @test (r₀ = [1., 0., 0.]; nout = [1., 2., 0.]; volfrac(box, nout, r₀) ≈ 0.25)
        @test (r₀ = [1., 0., 0.]; nout = [1., 2., 0.]; volfrac(box, -nout, r₀) ≈ 0.75)
    end

    @testset "quadrangular cylinder" begin
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

    @testset "general cases" begin
        # Test random cases.
        @test begin
            result = true
            for i = 1:100
                box = sort(randn(3,2), 2)
                r₀ = mean(box, 2)[:,1]
                nout = randn(3)
                result &= volfrac(box, nout, r₀) ≈ 0.5
            end
            result
        end

        @test begin
            result = true
            for i = 1:100
                box = sort(randn(3,2), 2)
                r₀ = randn(3)
                nout = randn(3)
                result &= (volfrac(box, nout, r₀) + volfrac(box, -nout, r₀) ≈ 1.)
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
