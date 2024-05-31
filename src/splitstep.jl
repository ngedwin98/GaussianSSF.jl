abstract type Integrator end

struct RK4IP <: Integrator
    dt::Float64
end

struct GSSFSim{M<:Integrator,T<:Model,S,C,N,D}
    method::M
    model::T
    N::Int
    L::Float64
    state::S
    cache::C
    linear_plans::D
    nonlinear_data::N
end

function GSSFSim(method::RK4IP, model::Model, N::Int, L::Real, Ω_funs)
    linear_plans = _init_plans(method, model, Ω_funs, N, L)
    state = zeros.(ComplexF64, size.(linear_plans))
    cache = (similar.(state), similar.(state), similar.(state))
    nonlinear_data = nonlinear_setup(model, L/N, method.dt, state)
    return GSSFSim(method, model, N, L, state, cache, linear_plans, nonlinear_data)
end

function _init_plans(method::RK4IP, model::Model, Ωk, N, L)
    k, _ = wavespace(N, L)
    k = ifftshift(k)
    cisdϕ = [cis.(-Ω.(2π.*k) .* method.dt ./ 2) for Ω ∈ Ωk]
    return linear_setup(model, cisdϕ...)
end

function init_state!(sim::GSSFSim, φ0_funs...)
    z, dz = realspace(sim.N, sim.L)
    for (i,φ0) ∈ enumerate(φ0_funs)
        x = sim.state[i]
        copyto!(x, [φ0(Z...)*sqrt(dz) for Z in product(repeated(z,ndims(x))...)])
    end
end

function init_state!(sim::GSSFSim{<:Integrator,<:ParallelMC}, φ0_funs...)
    z, dz = realspace(sim.N, sim.L)
    for (i,φ0) ∈ enumerate(φ0_funs)
        sim.state[i] .= Array(randn(ComplexF64, size(sim.state[i]))) ./ sqrt(2)
        sim.state[i] .+= φ0.(z) .* sqrt(dz)
    end
end

function step!(sim::GSSFSim{RK4IP}, t::Real)
    x = sim.state
    X, k, y = sim.cache
    nonlinear_step!(sim.model, k, x, t, sim.nonlinear_data...)
    linear_step!(sim.model, x, t, sim.linear_plans)
    copyto!.(X, x)
    linear_step!(sim.model, k, t, sim.linear_plans)
    _rk_substep!.(x, x, k, 6)
    _rk_substep!.(y, X, k, 2)
    nonlinear_step!(sim.model, k, y, t, sim.nonlinear_data...)
    _rk_substep!.(x, x, k, 3)
    _rk_substep!.(y, X, k, 2)
    nonlinear_step!(sim.model, k, y, t, sim.nonlinear_data...)
    _rk_substep!.(x, x, k, 3)
    _rk_substep!.(y, X, k, 1)
    linear_step!(sim.model, y, t, sim.linear_plans)
    nonlinear_step!(sim.model, k, y, t, sim.nonlinear_data...)
    linear_step!(sim.model, x, t, sim.linear_plans)
    _rk_substep!.(x, x, k, 6)
end

_rk_substep!(x, y, z, c) = (x .= y .+ z ./ c)


function gssf!(sim::GSSFSim{RK4IP}, Nsteps::Int, Nsave::Int=1, save_fun! = sim->Array.(sim.state))
    time = 0
    savesteps = Nsave == 1 ? [Nsteps+1] : round.(Int, range(1, Nsteps+1, length=Nsave))
    tout = (savesteps .-1).*sim.method.dt
    output = []
    
    next_save = popfirst!(savesteps)
    for t in 1 : Nsteps
        if t == next_save
            push!(output, save_fun!(sim))
            next_save = popfirst!(savesteps)
        end
        step!(sim, time)
        time += sim.method.dt
    end
    push!(output, save_fun!(sim))
    
    return output, tout
end

function wavespace(N, L)
    dξ = 1 / L
    ξ = fftfreq(N, N * dξ)
    return fftshift(ξ), dξ
end

function realspace(N, L)
    dz = L / N
    z = (-N/2 : N/2-1) * dz
    return collect(z), dz
end
