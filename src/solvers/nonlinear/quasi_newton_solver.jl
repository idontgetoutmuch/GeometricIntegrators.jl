
immutable QuasiNewtonSolver{T, TF, TJ, TL} <: AbstractNewtonSolver{T}
    @newton_solver_variables
    function QuasiNewtonSolver(x, Fparams, Jparams, linear_solver, nmax, atol, rtol, stol)
        new(x, Fparams, Jparams, linear_solver, NonlinearSolverParameters{T}(nmax, atol, rtol, stol), NonlinearSolverStatus{T}())
    end
end


function QuasiNewtonSolver(x::Vector, Fparams::NonlinearFunctionParameters; J=nothing, linear_solver=nothing, nmax=DEFAULT_nmax, atol=DEFAULT_atol, rtol=DEFAULT_rtol, stol=DEFAULT_stol, ϵ=DEFAULT_ϵ, autodiff=false)
    T = eltype(x)
    n = length(x)
    Jparams = getJacobianParameters(J, Fparams, ϵ, T, n, autodiff)
    linear_solver = getLinearSolver(T, n, linear_solver)
    QuasiNewtonSolver{T, typeof(Fparams), typeof(Jparams), typeof(linear_solver)}(x, Fparams, Jparams, linear_solver, nmax, atol, rtol, stol)
end


function solve!{T}(s::QuasiNewtonSolver{T})
    function_stages!(s.x, s.linear.b, s.Fparams)
    scale!(s.linear.b, -one(T))
    s.status.i  = 0
    s.status.rₐ = residual_absolute(s.linear.b)
    s.status.r₀ = s.status.rₐ

    if s.status.r₀ ≥ s.params.atol²
        computeJacobian(s.x, s.linear.A, s.Jparams)
        factorize!(s.linear)
        for s.status.i = 1:s.params.nmax
            solve!(s.linear)
            simd_xpy!(s.linear.b, s.x)
            s.status.rᵣ = residual_relative(s.linear.b, s.x)
            function_stages!(s.x, s.linear.b, s.Fparams)
            scale!(s.linear.b, -one(T))
            s.status.rₛ = s.status.rₐ
            s.status.rₐ = residual_absolute(s.linear.b)
            s.status.rₛ = abs(s.status.rₛ - s.status.rₐ)/s.status.r₀

            if s.status.rₐ < s.params.atol² || s.status.rᵣ < s.params.rtol || s.status.rₛ < s.params.stol²
                break
            end
        end
    end
    nothing
end
