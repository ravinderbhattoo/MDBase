export AbstractSimObj, MDSim, simulate, integrator, _stage, sim_init, MDParams, SParams, MParams, AbstractMDParams, IntgParams

export get_acceleration, set_acceleration!, get_potential_energy

export @SParams_generic_fields, @MParams_generic_fields

abstract type AbstractSimObj end
abstract type AbstractMDParams end

function Base.show(stream::IO, p::T) where T <: Union{AbstractSimObj,AbstractMDParams}
    println(stream, "Abstract Simulation/MD Object:")
    names = fieldnames(typeof(p))
    for n in names
        print(stream, "\t"), show(stream, n); println(stream)
    end
end

@def SParams_generic_fields begin
    N::I
    sim::SimObj
    kb::F
    mf_acc::F
    acc::Array{F, 3}
end

struct SParams{I <: IntType, SimObj <: AbstractSimObj, F <: FloatType} <: AbstractMDParams
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

function get_acceleration(v::Array{T,2}, u::Array{T,2}, params) where {T <: FloatType}
    dv = zeros(Types.F, size(v))
    set_acceleration!(dv, v, u, params)
    dv
end

function set_acceleration!(dv::Array{T,2}, v::Array{T,2}, u::Array{T,2}, params) where {T <: FloatType}
    for pot in params.S.sim.interatomic_potentials
        acceleration!(dv, v, u, pot, params)
    end
end

function get_potential_energy(v::Array{T,2}, u::Array{T,2}, params) where {T <: FloatType}
    pe = 0.0
    for pot in params.S.sim.interatomic_potentials
        pe += potential_energy(v, u, pot, params)
    end
    pe
end

struct MDSim{T1 <: FloatType, T2 <: IntType, T3<:Any, T4 <: SimulationBoundaries, pType <: AbstractMDParams} <: AbstractSimObj
    u0::Array{T1, 2}
    v0::Array{T1, 2}
    mass::Array{T1, 1}
    a_ids::Array{T2, 1}
    m_ids::Array{T2, 1}
    Δτ::T1
    save_every::T2
    thermo_save_every::T2
    interatomic_potentials::T3
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

    interatomic_potentials_ = MDBase.collect_interatomic_potentials(interatomic_potentials)

    MDSim(u0, v0, mass, a_ids, m_ids, Types.F(Δτ), save_every, thermo_save_every, interatomic_potentials_, boundary_condition, others)
end

function _stage(n::Types.I, sim::T, ensemble::Array{T2,1}; soode=nothing, verbose::Bool=false) where {T <: AbstractSimObj, T2 <: Ensemble}
    tspan = (0.0*sim.Δτ, n*sim.Δτ)
    params = exe_at_start(n, sim)
    print_thermo_at_start(params, sim, verbose)

    function SOODE(dv, v, u, p, t)
        fill!(dv, 0.0)
        set_acceleration!(dv, v, u, params)
        for ens in ensemble
            ddu!(dv, v, u, params, t, ens)
        end
    end
    if soode==nothing
        f = SOODE
    else
        f = soode
    end
    p = Float64[]
    prob = SecondOrderODEProblem(f, sim.v0, sim.u0, tspan, p, callback=mdcallbackset(params, ensemble))
    return prob, sim.Δτ, tspan[1]:sim.save_every*sim.Δτ:tspan[2], params
end

function _simulate(prob::T, dt::F, saveat, params) where {T,F}
    sol = solve(prob, VelocityVerlet(), dt=dt, saveat=saveat)
    params.M.step = 0
    sol
end

function simulate(n::Types.I, sim::T, ensemble::Array{T2,1}; soode=nothing, verbose::Bool=false) where {T <: AbstractSimObj, T2 <: Ensemble}
    prob, dt, saveat, params = _stage(n, sim, ensemble, soode=soode, verbose=verbose)
    sol = solve(prob, VelocityVerlet(), dt=dt, saveat=saveat)
    params.M.step = 0
    sol, params
end

function integrator(n::Types.I, sim::T, ensemble::Array{T2,1}; verbose::Bool=false) where {T <: AbstractSimObj, T2 <: Ensemble}
    prob, dt, saveat, params = _stage(n, sim, ensemble, verbose=verbose)
    init(prob, VelocityVerlet(), dt=dt, cb=mdcallbackset()), params
end

function exe_at_start(n::Types.I, sim::T) where T <: AbstractSimObj
    if MDBase.States.FreeRun
        _exe_at_start(n, sim)
    else
        throw("exe_at_start not implemented for object :: $(typeof(sim)).")
    end
end

function _exe_at_start(n::Types.I, sim::T) where T <: AbstractSimObj
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
    S = SParams(n, sim, kb, mf_acc, acc)
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
