module MDBase

export FloatType, IntType

using OrdinaryDiffEq: solve, init, step!, VelocityVerlet, DiscreteCallback, CallbackSet, SecondOrderODEProblem, get_du, ODESolution
using Reexport

@reexport using StaticArrays

abstract type CallBackType end
abstract type EoMType end

FloatType = Union{Float64, Float32}
IntType = Union{Int64}

mutable struct _States
    FreeRun::Bool
end

States = _States(false)

include("customtypes.jl")
include("boundaries.jl")
include("callbacks.jl")
include("potentials.jl")
include("ensembles.jl")
include("simulator.jl")
include("neighbors.jl")
include("io.jl")

end # module
