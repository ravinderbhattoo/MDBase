export AbstractSimObj, MDSim, simulate, integrator, sim_init, MDParams, SParams, MParams, AbstractMDParams, IntgParams

export get_acceleration, get_potential_energy

export @SParams_generic_fields, @MParams_generic_fields

abstract type AbstractSimObj end
abstract type AbstractMDParams end

@def SParams_generic_fields begin
    sim::SimObj
    N::I
    kb::F
    mf_acc::F
    acc::Array{F, 3}
end

struct SParams{SimObj <: AbstractSimObj, I <: IntType, F <: FloatType} <: AbstractMDParams
    @SParams_generic_fields
end

@def MParams_generic_fields begin
    step::I
    Temperature::F
    ke::F
end

mutable struct MParams{I <: IntType, F <: FloatType} <: AbstractMDParams
    @MParams_generic_fields
end

struct MDParams{T1 <: AbstractMDParams, T2 <: AbstractMDParams} <: AbstractMDParams
    S::T1
    M::T2
end

struct NullMDParams <: AbstractMDParams
end

struct IntgParams{T1 <: Ensemble, T2 <: AbstractMDParams} <: AbstractMDParams
    ensemble::Array{T1, 1}
    params::T2
end

function SOODE(v, u, p, t)
    ensemble, params = p[1].ensemble, p[1].params
    dv = get_acceleration(v, u, params, params.S.sim)
    for ens in ensemble
        ddu!(dv, v, u, params, t, ens)
    end
    return dv
end

function get_acceleration(v::Array{T,2}, u::Array{T,2}, params, sim::SimObj) where {T <: FloatType, SimObj <: AbstractSimObj}
    dv = zeros(Types.F, size(v))
    for pot in params.S.sim.interatomic_potentials
        acceleration!(dv, v, u, pot, params)
    end
    dv
end

function get_potential_energy(u::Array{T,2}, params, sim::SimObj) where {T <: FloatType, SimObj <: AbstractSimObj}
    pe = 0.0
    for pot in params.S.sim.interatomic_potentials
        pe += potential_energy(u, pot, params)
    end
    pe
end

struct MDSim{T1 <: FloatType, T2 <: IntType, T3 <: PotentialParameters, T4 <: SimulationBoundaries, pType <: AbstractMDParams} <: AbstractSimObj
    u0::Array{T1, 2}
    v0::Array{T1, 2}
    mass::Array{T1, 1}
    a_ids::Array{T2, 1}
    m_ids::Array{T2, 1}
    Δτ::T1
    save_every::T2
    thermo_save_every::T2
    interatomic_potentials::Array{T3, 1}
    boundary_condition::T4
    others::pType
end

function MDSim(u0, v0, mass, interatomic_potentials, boundary_condition; a_ids = nothing, m_ids = nothing, Δτ = Types.F(1.0), save_every = Types.I(100), thermo_save_every = Types.I(100), others=NullMDParams())

    if a_ids==nothing
        a_ids = ones(Types.I, size(u0,2))
        m_ids = copy(a_ids)
    elseif m_ids==nothing
        m_ids = copy(a_ids)
    end

    MDSim(u0, v0, mass, a_ids, m_ids, Types.F(Δτ), save_every, thermo_save_every, interatomic_potentials, boundary_condition, others)
end


function simulate(Sim::T, n::Types.I, ensemble::Array{T2,1}; verbose::Bool=false) where {T <: AbstractSimObj, T2 <: Ensemble}
    tspan = (0.0*Sim.Δτ, n*Sim.Δτ)
    params = exe_at_start(Sim, n)
    print_thermo_at_start(params, Sim, verbose)
    p = IntgParams(ensemble, params)
    prob = SecondOrderODEProblem(SOODE, Sim.v0, Sim.u0, tspan, p, callback=mdcallbackset())
    solve(prob, VelocityVerlet(), dt=Sim.Δτ, saveat=tspan[1]:Sim.save_every*Sim.Δτ:tspan[2])
end

function integrator(Sim::T, n::Types.I, ensemble::Array{T2,1}; verbose::Bool=false) where {T <: AbstractSimObj, T2 <: Ensemble}
    tspan = (0.0*Sim.Δτ, n*Sim.Δτ)
    params = exe_at_start(Sim, n)
    print_thermo_at_start(params, Sim, verbose)
    p = IntgParams(ensemble, params)
    prob = SecondOrderODEProblem(SOODE, Sim.v0, Sim.u0, tspan, p, callback=mdcallbackset())
    init(prob, VelocityVerlet(), dt=Sim.Δτ)
end

function exe_at_start(sim::T, n::Types.I) where T <: AbstractSimObj
    if MDBase.States.FreeRun
        _exe_at_start(sim, n)
    else
        throw("exe_at_start not implemented for object :: $(typeof(sim)).")
    end
end

function _exe_at_start(sim::T, n::Types.I) where T <: AbstractSimObj
    ux = @view sim.u0[1, :]
    uy = @view sim.u0[2, :]
    uz = @view sim.u0[3, :]
    vx = @view sim.v0[1, :]
    vy = @view sim.v0[2, :]
    vz = @view sim.v0[3, :]
    apply_simulation_bc!(ux, uy, uz, vx, vy, vz, sim.boundary_condition)
    N = Types.I(size(sim.u0, 2))
    kb = Types.F(1.0)
    mf_acc = Types.F(1.0)
    acc = zeros(Types.F, (size(sim.u0)..., n+1))
    S = SParams(sim, N, kb, mf_acc, acc)
    M = MParams(Types.I(0), Types.F(0), Types.F(0))
    return MDParams(S, M)
end

function print_thermo_at_start(params, sim::AbstractSimObj, verbose::Bool)
    if MDBase.States.FreeRun
        @warn "Thermodynamic print is empty. No thermo values to print. Please write a print_thermo(params, sim::AbstractSimObj, verbose::Bool) function."
    else
        throw("print_thermo not implemented for object :: $(typeof(sim)).")
    end
end


#
# function sim_init(sim::AbstractSimObj; kwargs...)
#     args = []
#     k = keys(kwargs)
#     for f in fieldnames(typeof(sim))
#         if f in k
#             push!(args, kwargs[Symbol(f)])
#         else
#             push!(args, getproperty(sim,f))
#         end
#     end
#     AbstractSimObj(args...)
# end
#
# function sim_init(res::DiffEqBase.ODESolution, sim::AbstractSimObj)
#     data = Array(reshape(Array(res.u[end]),(3,:)))
#     N = Int(size(data, 2)/2)
#     v0 = data[:, 1:N]
#     u0 = data[:, N+1:end]
#     sim_init(sim, u0=u0, v0=v0)
# end
