struct LinearPlan{T,P1,P2}
    forward!::P1
    scalars::T
    backward!::P2
end

Base.size(plan::LinearPlan) = size(plan.scalars)

# Perform 1D linear operation for 1D fields
LinearPlan(dD::AbstractArray) = LinearPlan(plan_fft!(dD), dD/length(dD), plan_bfft!(dD))

# Perform 2D linear operation for symmetric fields
LinearPlan(dD1::AbstractVector, dD2::AbstractVector) = LinearPlan(dD1 * transpose(dD2))

# Perform 2D linear operation via 2x 1D FFTs for nonsymmetric Hermitian fields
function LinearPlan(dD1::Adjoint, dD2::AbstractVector)
    # dD = transpose(dD1) * transpose(dD2)
    dD = conj.(dD1.parent) * transpose(dD2)
    forward! = (plan_bfft!(dD,1), plan_fft!(dD,2))
    backward! = (plan_fft!(dD,1), plan_bfft!(dD,2))
    return LinearPlan(forward!, dD/length(dD), backward!)
end

# Tuple construct for 2D FFT as two 1D FFTs (for nonsymmetric Hermitian fields)
function Base.:(*)(ft_plan!::Tuple{AbstractFFTs.Plan,AbstractFFTs.Plan}, x)
    for plan! in ft_plan!
        plan! * x
    end
end

function LinearPlan(dD::AbstractVector, Ntraj::Int)
    N = length(dD)
    dD = hcat((copy(dD) for i in 1:Ntraj)...)
    return LinearPlan(plan_fft!(dD,1), dD/N, plan_bfft!(dD,1))
end

function taylor(Ω_coeffs::Vector{<:Real}=[0,0,1])
    return let Ω_coeffs=Ω_coeffs
        k -> sum(Ωi/factorial(i-1)*k.^(i-1) for (i,Ωi) in enumerate(Ω_coeffs))
    end
end
taylor(Ω_coeffs...) = taylor(collect(Ω_coeffs))
