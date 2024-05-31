abstract type Model end

function linear_setup(::Model, cisdϕ...)
    dD = LinearPlan.(cisdϕ)
    for i ∈ eachindex(cisdϕ)
        for j ∈ i : length(cisdϕ)
            dD = (dD...,
                LinearPlan(cisdϕ[i],cisdϕ[j]), LinearPlan(adjoint(cisdϕ[i]),cisdϕ[j])
            )
        end
    end
    return dD
end

function linear_step!(::Model, state, t::Real, plans)
    for (x, plan) ∈ zip(state, plans)
        plan.forward! * x
        x .*= plan.scalars
        plan.backward! * x
    end
end

abstract type DerivedModel <: Model end
nonlinear_setup(derived::DerivedModel, args...) = nonlinear_setup(derived.model, args...)
linear_setup(derived::DerivedModel, args...) = linear_setup(derived.model, args...)

struct QD3WM <: Model
    ϵ::Float64
end
nonlinear_setup(model::QD3WM, dz, dt, x0) = (model.ϵ / sqrt(dz) * dt,)
function nonlinear_step!(::QD3WM, x′, x, t, dϵt)
    a, b, aa, a⁺a, ab, a⁺b, bb, b⁺b = x
    a′, b′, aa′, a⁺a′, ab′, a⁺b′, bb′, b⁺b′ = x′
      a′ .= b .* conj.(a) .+ diag(a⁺b)
      b′ .= -1/2 .* (a .* a .+ diag(aa))
     aa′ .= b .* a⁺a .+ conj.(a) .* transpose(ab) .+ transpose(b) .* (conj.(a⁺a) + I) .+ adjoint(a) .* ab
    a⁺a′ .= conj.(b) .* aa .+ a .* adjoint(a⁺b) .+ transpose(b) .* conj.(aa) .+ adjoint(a) .* a⁺b
     ab′ .= b .* a⁺b .+ conj.(a) .* bb .- transpose(a) .* aa
    a⁺b′ .= conj.(b) .* ab .+ a .* b⁺b .- transpose(a) .* a⁺a
     bb′ .= .-(a .* ab .+ transpose(a) .* transpose(ab))
    b⁺b′ .= .-(conj.(a) .* a⁺b .+ transpose(a) .* adjoint(a⁺b))
    apply_scalar!(x′, dϵt)
end

struct NLSE <: Model
    g::Float64
end
nonlinear_setup(model::NLSE, dz, dt, x0) = (1im*model.g / dz * dt, similar(x0[1]))
function nonlinear_step!(::NLSE, x′, x, t, idgt, u)
    a, aa, a⁺a = x
    a′, aa′, a⁺a′ = x′
    a′ .= conj.(a) .* diag(aa)
    u .= abs2.(a) .+ diag(a⁺a)
    a′ .+= a .* (u .+ diag(a⁺a))
    aa′ .= 2 .* (u .+ transpose(u)) .* aa
    a⁺a′ .= -2 .* (u .- transpose(u)) .* a⁺a
    u .= a.*a .+ diag(aa)
    aa′ .+= transpose(u) .* (conj.(a⁺a) + I) .+ u .* a⁺a
    a⁺a′ .+= transpose(u) .* conj.(aa) .- conj.(u) .* aa
    apply_scalar!(x′, -idgt)
end

struct Squeezing{M<:Model} <: DerivedModel
    model::M
end

function linear_setup(::Squeezing{QD3WM}, cisdϕ_a, cisdϕ_b)
    LinearPlan(cisdϕ_b), LinearPlan(cisdϕ_a,cisdϕ_a),
        LinearPlan(adjoint(cisdϕ_a),cisdϕ_a)
end
function nonlinear_step!(::Squeezing{QD3WM}, x′, x, t, dϵt)
    b, aa, a⁺a = x
    b′, aa′, a⁺a′ = x′
    b′ .= -1/2 .* diag(aa)
    aa′ .= b .* a⁺a .+ transpose(b) .* (conj.(a⁺a) + I)
    a⁺a′ .= conj.(b) .* aa .+ transpose(b) .* conj.(aa)
    apply_scalar!(x′, dϵt)
end

struct Linearized{M<:Model} <: DerivedModel
    model::M
end

function nonlinear_step!(::Linearized{Squeezing{QD3WM}}, x′, x, t, dϵt)
    b, aa, a⁺a = x
    b′, aa′, a⁺a′ = x′
    b′ .= 0
    aa′ .= b .* a⁺a .+ transpose(b) .* (conj.(a⁺a) + I)
    a⁺a′ .= conj.(b) .* aa .+ transpose(b) .* conj.(aa)
    apply_scalar!(x′, dϵt)
end

function nonlinear_step!(::Linearized{NLSE}, x′, x, t, idgt, u)
    a, aa, a⁺a = x
    a′, aa′, a⁺a′ = x′
    a′ .= a .* (abs2.(a))
    u .= abs2.(a)
    aa′ .= 2 .* (u .+ transpose(u)) .* aa
    a⁺a′ .= -2 .* (u .- transpose(u)) .* a⁺a
    u .= a.*a
    aa′ .+= transpose(u) .* (conj.(a⁺a) + I) .+ u .* a⁺a
    a⁺a′ .+= transpose(u) .* conj.(aa) .- conj.(u) .* aa
    apply_scalar!(x′, -idgt)
end

struct Classical{M<:Model} <: DerivedModel
    model::M
end
linear_setup(::Classical, cisdϕ...) = LinearPlan.(cisdϕ)
function nonlinear_step!(::Classical{QD3WM}, x′, x, t, dϵt)
    a, b = x
    a′, b′ = x′
    a′ .= b .* conj.(a)
    b′ .= -1/2 .* a.*a
    apply_scalar!(x′, dϵt)
end
function nonlinear_step!(::Classical{NLSE}, x′, x, t, idgt, u)
    a, = x
    a′, = x′
    a′ .= a .* abs2.(a)
    apply_scalar!(x′, -idgt)
end

struct ParallelMC{M<:Classical} <: DerivedModel
    model::M
    Ntraj::Int
end
linear_setup(par::ParallelMC, cisdϕ...) where T = LinearPlan.(cisdϕ, par.Ntraj)
nonlinear_step!(par::ParallelMC, data...) where T = nonlinear_step!(par.model, data...)

function apply_scalar!(x, h)
    for xi in x
        xi .*= h
    end
end
