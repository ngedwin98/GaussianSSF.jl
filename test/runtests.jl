using GaussianSSF
using Test
using FFTW          # concrete CPU FFT backend for `Array`; required for any test to run
using LinearAlgebra # norm / diag / Diagonal / I for the CW gain-spectrum benchmark

# ---------------------------------------------------------------------------
# Backends under test.
#
# The CPU backend (`Array`) is always exercised. The GPU backend (`CuArray`) is
# added only when CUDA is installed *and* a functional device is present, so the
# default `Pkg.test()` on a CPU-only machine pulls in no CUDA stack. CUDA is
# deliberately not a declared dependency of the package or its tests.
# ---------------------------------------------------------------------------
const BACKENDS = Any[Array]

if Base.find_package("CUDA") !== nothing
    try
        @eval using CUDA
        if CUDA.functional()
            push!(BACKENDS, CuArray)
        else
            @info "CUDA is installed but not functional; skipping GPU tests."
        end
    catch err
        @info "CUDA present but failed to load; skipping GPU tests." exception = err
    end
end

# True only when a functional GPU backend was added above. Used instead of the
# `CuArray` symbol, which is undefined when CUDA was never loaded.
const HAS_GPU = length(BACKENDS) > 1

@info "Testing GaussianSSF on backends: $(BACKENDS)"

# Bring a (possibly device-resident) state array back to a plain host Array.
host(x) = Array(x)

@testset "GaussianSSF" begin

    @testset "Grid constructors" begin
        N, L = 32, 10.0
        x, dx = realspace(N, L)
        ξ, dξ = wavespace(N, L)
        @test length(x) == N
        @test dx ≈ L / N
        @test dξ ≈ 1 / L
        @test x[1] ≈ -L / 2                # grid starts at -L/2
        @test x[end] ≈ L / 2 - dx          # ...and excludes +L/2 (aliased with x[1] under periodicity)
        @test issorted(x)
        @test step(x[1]:dx:x[end]) ≈ dx    # uniform spacing
        @test length(ξ) == N
        @test dx * dξ ≈ 1 / N              # exact DFT reciprocity: real & reciprocal grids are conjugate
    end

    # -----------------------------------------------------------------------
    # Real and reciprocal grids form an exact FFT conjugate pair.
    #
    # This is the invariant behind the two historical grid fixes (ba70914,
    # dbed4cd): a plane wave exp(2πi m x/L) sampled on the real grid must
    # transform to a *single* reciprocal-grid mode at ξ = m/L, with no leakage.
    # It holds iff dx·dξ = 1/N exactly and the ifftshift/fftshift centering
    # agrees — the same convention `_init_plans` relies on to build the
    # dispersion propagators. If the endpoint-inclusive spacing (L/(N-1)) ever
    # crept back in, the mode would smear across bins and this would fail.
    # -----------------------------------------------------------------------
    @testset "Grids form an FFT conjugate pair" begin
        for (N, L) in ((32, 10.0), (64, 8.0))
            x, dx = realspace(N, L)
            ξ, _ = wavespace(N, L)
            for m in (1, 3, -5)
                f = cis.(2π * m .* x ./ L)
                F = fftshift(fft(ifftshift(f)))
                peak = argmax(abs.(F))
                @test ξ[peak] ≈ m / L                # energy lands on the correct mode
                @test abs(F[peak]) ≈ sum(abs, F)     # ...and *only* that mode (no leakage)
            end
        end
    end

    @testset "Backend: $bk" for bk in BACKENDS

        # -------------------------------------------------------------------
        # FFT roundtrip + normalization.
        #
        # With zero dispersion (D ≡ 0) and zero nonlinearity, a single RK4IP
        # step must be the identity: forward FFT → (×1) → backward FFT with the
        # 1/N normalization must recover the input exactly. This validates that
        # the chosen FFT backend is loaded and that grid ordering / scaling are
        # consistent.
        # -------------------------------------------------------------------
        @testset "FFT roundtrip is identity (D=0, g=0)" begin
            N, L = 32, 8.0
            model = Classical(NLSE(0.0))              # nonlinearity switched off
            D = taylor(0.0)                           # dispersion switched off
            sim = GSSFSim(RK4IP(0.01), model, N, L, (D,); arrT = bk)
            x, _ = realspace(N, L)
            # arbitrary complex profile
            init_state!(sim, xi -> cospi(xi / L) + 0.3im * sinpi(2xi / L))
            initial = host(sim.state[1])
            step!(sim, 0.0)
            @test host(sim.state[1]) ≈ initial rtol = 1e-12
        end

        # -------------------------------------------------------------------
        # Pure dispersion conserves total photon number.
        #
        # With real dispersion and no nonlinearity, propagation is unitary in
        # reciprocal space, so Σ|ψ|² (the discrete photon number) is invariant.
        # -------------------------------------------------------------------
        @testset "Pure dispersion conserves photon number" begin
            N, L = 64, 10.0
            β2 = 2.0
            model = Classical(NLSE(0.0))
            D = taylor(0.0, 0.0, β2)                  # Ω(κ) = β2/2 κ²
            Z, Nsteps = 0.5, 200
            sim = GSSFSim(RK4IP(Z / Nsteps), model, N, L, (D,); arrT = bk)
            init_state!(sim, xi -> exp(-xi^2 / 2) * cis(0.7xi))
            n0 = sum(abs2, host(sim.state[1]))
            gssf!(sim, Nsteps)
            n1 = sum(abs2, host(sim.state[1]))
            @test n1 ≈ n0 rtol = 1e-10
        end

        # -------------------------------------------------------------------
        # Pure dispersion: exact phase on a single Fourier mode.
        #
        # A grid-aligned plane wave exp(2πi m x / L) is an eigenmode of the
        # linear step. After propagating a distance Z under Ω(κ) = β2/2 κ², it
        # acquires exactly the phase exp(-i Ω(2π m/L) Z) and is otherwise
        # unchanged. This pins down the dispersion sign, the 2π factor, and the
        # ifftshift/fftshift grid conventions.
        # -------------------------------------------------------------------
        @testset "Pure dispersion: exact single-mode phase" begin
            N, L = 64, 10.0
            β2 = 2.0
            m = 3                                     # integer mode index
            model = Classical(NLSE(0.0))
            D = taylor(0.0, 0.0, β2)
            Z, Nsteps = 0.5, 50
            sim = GSSFSim(RK4IP(Z / Nsteps), model, N, L, (D,); arrT = bk)
            init_state!(sim, xi -> cis(2π * m * xi / L))
            initial = host(sim.state[1])
            gssf!(sim, Nsteps)
            θ = (β2 / 2) * (2π * m / L)^2 * Z         # Ω(2π m/L) · Z
            expected = initial .* cis(-θ)
            @test host(sim.state[1]) ≈ expected rtol = 1e-10
        end

        # -------------------------------------------------------------------
        # Classical self-phase modulation (SPM).
        #
        # A spatially uniform classical field with no dispersion obeys
        # da/dz = -i g |a|² a, i.e. it rotates in phase at rate g|A|² with |a|
        # fixed: a(Z) = A exp(-i g |A|² Z).
        # -------------------------------------------------------------------
        @testset "Classical SPM phase rotation" begin
            N, L = 16, 4.0
            g, A = 1.0, 1.0
            model = Classical(NLSE(g))
            D = taylor(0.0)                           # no dispersion
            Z, Nsteps = 1.0, 2000
            sim = GSSFSim(RK4IP(Z / Nsteps), model, N, L, (D,); arrT = bk)
            init_state!(sim, _ -> complex(A))         # uniform field
            initial = host(sim.state[1])
            gssf!(sim, Nsteps)
            expected = initial .* cis(-g * A^2 * Z)
            @test host(sim.state[1]) ≈ expected rtol = 1e-4
        end

        # -------------------------------------------------------------------
        # Smoke tests: every model builds and runs with finite, well-shaped
        # output. Kept at small N so the 2D covariance FFTs stay cheap.
        # -------------------------------------------------------------------
        @testset "Smoke: models run and stay finite" begin
            N, L = 16, 8.0

            allfinite(state) = all(x -> all(isfinite, host(x)), state)

            # NLSE and its reductions (single field)
            for model in (NLSE(0.5), Linearized(NLSE(0.5)), Classical(NLSE(0.5)))
                sim = GSSFSim(RK4IP(0.005), model, N, L, (taylor(0, 0, 1.0),); arrT = bk)
                init_state!(sim, xi -> exp(-xi^2))
                out, _ = gssf!(sim, 10)
                @test allfinite(sim.state)
                @test allfinite(out[end])
            end

            # QD3WM (two mean fields: signal a + pump b). State = (a, b, Σ/Π…);
            # both means are vectors, so seed signal = vacuum, pump = a pulse.
            Dfh = taylor(0, 0, 1.0)
            Dsh = taylor(0, 0, 2.0)
            for model in (QD3WM(1.0), Classical(QD3WM(1.0)))
                sim = GSSFSim(RK4IP(0.005), model, N, L, (Dfh, Dsh); arrT = bk)
                init_state!(sim, zero, xi -> 5exp(-xi^2 / 2))
                out, _ = gssf!(sim, 10)
                @test allfinite(sim.state)
                @test allfinite(out[end])
            end

            # Squeezing reductions eliminate the signal mean and pump covariances:
            # state = (pump SH mean, Σ11, Π11). Only field #1 (the pump) is a
            # vector to seed; the signal covariances start from vacuum (zero).
            for model in (Squeezing(QD3WM(1.0)), Linearized(Squeezing(QD3WM(1.0))))
                sim = GSSFSim(RK4IP(0.005), model, N, L, (Dfh, Dsh); arrT = bk)
                init_state!(sim, xi -> 5exp(-xi^2 / 2))
                out, _ = gssf!(sim, 10)
                @test allfinite(sim.state)
                @test allfinite(out[end])
            end

            # Parallel Monte-Carlo classical trajectories
            mc = ParallelMC(Classical(NLSE(0.5)), 8)
            sim = GSSFSim(RK4IP(0.005), mc, N, L, (taylor(0, 0, 1.0),); arrT = bk)
            init_state!(sim, xi -> exp(-xi^2))
            @test size(sim.state[1]) == (N, 8)        # columns = trajectories
            out, _ = gssf!(sim, 10)
            @test allfinite(sim.state)
        end
    end

    # -----------------------------------------------------------------------
    # Element type (precision) selection via `elT`.
    #
    # Default is ComplexF64; ComplexF32 is accepted and propagates to the state
    # (and, internally, to the FFT plans, so a full run stays in F32). Only
    # complex eltypes are allowed — the constructor's `elT::Type{<:Complex}`
    # bound rejects real types at the call boundary. CPU-only; backend-agnostic.
    # -----------------------------------------------------------------------
    @testset "Element type (elT)" begin
        N, L = 16, 8.0
        model = Classical(NLSE(0.5))
        D = (taylor(0, 0, 1.0),)

        # default precision
        sim64 = GSSFSim(RK4IP(0.005), model, N, L, D)
        @test eltype(sim64.state[1]) == ComplexF64

        # single precision runs end-to-end and stays ComplexF32 throughout
        sim32 = GSSFSim(RK4IP(0.005), model, N, L, D; elT = ComplexF32)
        @test eltype(sim32.state[1]) == ComplexF32
        init_state!(sim32, xi -> exp(-xi^2))
        out, _ = gssf!(sim32, 10)
        @test all(x -> all(isfinite, host(x)), out[end])
        @test eltype(out[end][1]) == ComplexF32

        # real element types are rejected by the `elT::Type{<:Complex}` bound
        @test_throws TypeError GSSFSim(RK4IP(0.005), model, N, L, D; elT = Float64)
    end

    # -----------------------------------------------------------------------
    # CW (flat-pump) parametric squeezing — exact closed form, both sectors.
    # See docs/squeezing_benchmark.md.
    #
    # A translationally-invariant (flat) SH pump makes undepleted squeezing
    # diagonal in frequency: each Fourier mode ω pairs only with -ω to form an
    # independent two-mode squeezer with generator [D γ; -γ -D], D(ω)=½sω²,
    # γ=εβ₀.  From vacuum the Bogoliubov coefficients are
    #   A(ω) = cosh(νZ) - i(D/ν) sinh(νZ),   B(ω) = -i(γ/ν) sinh(νZ),  ν=√(γ²-D²)
    # (complex ν ⇒ amplification for D<γ and sin/cos sidelobes for D>γ), giving
    #   n(ω) = ⟨a†_ω a_ω⟩ = |B|²          (photon-number / gain spectrum)
    #   m(ω) = ⟨a_ω a_{-ω}⟩ = i·A·B        (entangled-pair correlation; i = pkg phase)
    # and n(0)=sinh²(γZ).  On the discrete grid each mode is *exactly* such a
    # squeezer, so in the Fourier basis Π=⟨a†a⟩ is DIAGONAL (each mode marginally
    # thermal) while Σ=⟨aa⟩ is ANTI-DIAGONAL (ω correlated only with -ω).  Both
    # match the closed form to the RK4IP floor (~1e-8) — a strong, non-monotonic,
    # entanglement-sensitive check.  CPU-only; no GPU required.
    # -----------------------------------------------------------------------
    @testset "CW parametric squeezing (flat pump)" begin
        N, Lwin, s, ε, β0, Z, Nsteps = 128, 20.0, 1.0, 1.0, 1.0, 0.5, 500
        model = Linearized(Squeezing(QD3WM(ε)))
        sim = GSSFSim(RK4IP(Z / Nsteps), model, N, Lwin, (taylor(0.0, 0.0, s), taylor(0.0)))
        init_state!(sim, _ -> β0)               # flat pump ⇒ translational invariance
        gssf!(sim, Nsteps)
        Π = host(sim.state[3])
        Σ = host(sim.state[2])

        # Fourier basis: number sector ⟨a†_ω a_ω'⟩ = WΠW†; pair sector ⟨a_ω a_ω'⟩ = WΣWᵀ.
        W  = fft(Matrix{ComplexF64}(I, N, N), 1) / sqrt(N)   # unitary DFT
        Πf = W * Π * W'
        Σf = W * Σ * transpose(W)
        revω = [mod(-(i - 1), N) + 1 for i in 1:N]           # index of -ω for each ω
        n_num = real.(diag(Πf))                              # ⟨a†_ω a_ω⟩   = n(ω)
        m_num = [Σf[i, revω[i]] for i in 1:N]                # ⟨a_ω a_{-ω}⟩ = m(ω)

        ω = 2π .* fftfreq(N, N / Lwin)                       # native order, matches W
        D = (s / 2) .* ω .^ 2
        γ = ε * β0
        function theory(D)                                   # → (n(ω), m(ω))
            ν  = sqrt(complex(γ^2 - D^2))
            sc = iszero(ν) ? complex(Z) : sinh(ν * Z) / ν
            B  = -im * γ * sc
            A  = cosh(ν * Z) - im * D * sc
            return abs2(B), im * A * B
        end
        n_th = first.(theory.(D))
        m_th = last.(theory.(D))

        # Number sector: Π diagonal, gain spectrum n(ω) incl. sinh/sin sidelobes.
        @test norm(Πf - Diagonal(diag(Πf))) / norm(Πf) < 1e-10          # translational invariance
        @test n_num[1] ≈ sinh(γ * Z)^2 rtol = 1e-5                      # DC amplification anchor
        @test maximum(abs.(n_num .- n_th)) / maximum(n_th) < 1e-5       # full gain spectrum

        # Pair sector: Σ anti-diagonal (entanglement is ω ↔ -ω only), exact m(ω).
        Σanti = [i == revω[j] ? Σf[i, j] : zero(eltype(Σf)) for i in 1:N, j in 1:N]
        @test norm(Σf - Σanti) / norm(Σf) < 1e-10                       # only ω ↔ -ω correlations
        @test maximum(abs.(m_num .- m_th)) / maximum(abs.(m_th)) < 1e-5 # exact complex pair correlation
    end

    # -----------------------------------------------------------------------
    # CPU/GPU equivalence: identical code, identical results to tight tolerance.
    # Only runs when a functional CUDA device pushed CuArray into BACKENDS.
    # -----------------------------------------------------------------------
    @testset "CPU/GPU backend equivalence" begin
        if !HAS_GPU
            @test_skip true    # no functional GPU in this environment
        else
            N, L = 32, 10.0
            model = QD3WM(1.0)
            D = (taylor(0, 0, 1.0), taylor(0, 0, 2.0))
            Z, Nsteps = 0.2, 40
            results = map(BACKENDS) do bk
                sim = GSSFSim(RK4IP(Z / Nsteps), model, N, L, D; arrT = bk)
                init_state!(sim, zero, xi -> 5exp(-xi^2 / 2))
                gssf!(sim, Nsteps)[1][end]
            end
            cpu, gpu = results[1], results[end]
            for (c, g) in zip(cpu, gpu)
                @test Array(c) ≈ Array(g) rtol = 1e-8
            end
        end
    end
end
