export PotentialParameters, PairPotential, MLPotential
export Dummy
export acceleration!, potential_energy, potential_force

abstract type PotentialParameters end

abstract type MLPotential <: PotentialParameters end
abstract type PairPotential <: PotentialParameters end
abstract type _3BodyPotential <: PotentialParameters end
abstract type _4BodyPotential <: PotentialParameters end
abstract type NBodyPotential <: PotentialParameters end
abstract type BondedPotential <: PotentialParameters end
abstract type FieldPotentials <: PotentialParameters end

struct Dummy <: PotentialParameters
end

@inline function potential_energy(r, pot::T) where T <: PotentialParameters
    if MDBase.States.FreeRun
        0.0
    else
        throw("Acceleration is not implemented for potential :: $(typeof(pot)). Signature: acceleration!(dv, v, u, pot::$(typeof(pot)), params).")
    end
end

@inline function potential_force(r, pot::T) where T <: PotentialParameters
    if MDBase.States.FreeRun
        0.0
    else
        throw("Acceleration is not implemented for potential :: $(typeof(pot)). Signature: acceleration!(dv, v, u, pot::$(typeof(pot)), params).")
    end
end

function acceleration!(dv, v, u, pot::T, params) where T <: PotentialParameters
    if MDBase.States.FreeRun
        dv
    else
        throw("Acceleration is not implemented for potential :: $(typeof(pot)). Signature: acceleration!(dv, v, u, pot::$(typeof(pot)), params).")
    end
end

function potential_energy(v, u, pot::T, params) where T <: PotentialParameters
    if MDBase.States.FreeRun
        0.0
    else
        throw("Potential energy is not implemented for potential :: $(typeof(pot)). Signature: potential_energy(u, pot::$(typeof(pot)), params)")
    end
end
