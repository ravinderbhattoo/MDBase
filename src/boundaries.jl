export GeneralBC, BCTypes, PeriodicBoundary, ReflectiveBoundary, OpenBoundary
export CubicPBC, apply_simulation_bc!, applybc!, distance, lattice

"""
    SimulationBoundaries

Abstract type to define general functuonality for simulations boudaries.
"""
abstract type SimulationBoundaries end

"""
    BCTypes

Abstract type to define general functuonality for boundary condition types.
"""
abstract type BCTypes end

"""
    PeriodicBoundary{T} <: BCTypes where T

Define a periodic boundary.

# Filednames:
- `X0 :: xlow value`
- `X1 :: xhigh value`
- `L :: Length i.e X1 - X0`

# Examples
```jldoctest
julia> PeriodicBoundary(X0, X1, L)
```
"""
struct PeriodicBoundary{T} <: BCTypes where T
    X0::T
    X1::T
    L::T
end

function Base.show(stream::IO, BC::PeriodicBoundary)
    println(stream, "Periodic Boundary (Low: $(BC.X0) High: $(BC.X1) L: $(BC.L))")
end

"""
    PeriodicBoundary{T} <: BCTypes where T

Define a periodic boundary.

# Input Args:
- `X0 :: xlow value`
- `X1 :: xhigh value`

# Output:
- struct of type PeriodicBoundary

# Examples
```jldoctest
julia> PeriodicBoundary(X0, X1)
```
"""
function PeriodicBoundary(X0::T1, X1::T2) where {T1, T2}
    L = X1 - X0
    return PeriodicBoundary(promote(X0, X1, L)...)
end

"""
    ReflectiveBoundary{T} <: BCTypes where T

Define a reflective boundary.

# Filednames:
- `X0 :: xlow value`
- `X1 :: xhigh value`
- `L :: Length i.e X1 - X0`
- `skin :: to reduce penetration in reflective boundary`

# Examples
```jldoctest
julia> PeriodicBoundary(X0, X1, L)
```
"""
struct ReflectiveBoundary{T} <: BCTypes where T
    X0::T
    X1::T
    L::T
    skin::T
end

"""
    ReflectiveBoundary{T} <: BCTypes where T

Define a reflective boundary.

# Input Args:
- `X0 :: xlow value`
- `X1 :: xhigh value`
- `skin :: to reduce penetration in reflective boundary. (default 0.001)`

# Output:
- struct of type ReflectiveBoundary

# Examples
```jldoctest
julia> ReflectiveBoundary(X0, X1)
```
"""
function ReflectiveBoundary(X0::T1, X1::T2; skin=0.001) where {T1, T2}
    L = X1 - X0
    X0 += skin*L
    X1 -= skin*L
    L = X1 - X0
    return ReflectiveBoundary(promote(X0, X1, L, skin)...)
end

"""
    OpenBoundary

It doesn't require anything to impose constriants.
"""
struct OpenBoundary <: BCTypes
end

struct GeneralBC{T1 <: BCTypes, T2 <: BCTypes, T3 <: BCTypes} <: SimulationBoundaries
    X::T1
    Y::T2
    Z::T3
end

function Base.show(stream::IO, BC::GeneralBC)
    println(stream, "General(Orthogonal) Boundary Condition:")
    print(stream, "\tDirection 1: "); show(stream, BC.X)
    print(stream, "\tDirection 2: "); show(stream, BC.Y)
    print(stream, "\tDirection 3: "); show(stream, BC.Z)
end


function CubicPBC(L::T) where T
    X = PeriodicBoundary(0L,L,L)
    GeneralBC(X,X,X)
end

function apply_simulation_bc!(ux, uy, uz, vx, vy, vz, BC::GeneralBC)
    applybc!(ux, vx, BC.X)
    applybc!(uy, vy, BC.Y)
    applybc!(uz, vz, BC.Z)
end

function applybc!(x, v, BC::OpenBoundary)
    nothing
end

function applybc!(x, v, BC::PeriodicBoundary)
    # @inbounds for i in eachindex(x)
    #     if x[i]<BC.X0 || x[i]>BC.X1
    #         x[i] = BC.X0 - fld(x[i]-BC.X0, BC.L)*BC.L
    #     else
    #     end
    # end
    @. x = x - fld(x - BC.X0, BC.L) * BC.L
end

function applybc!(x, v, BC::ReflectiveBoundary)
    @inbounds for i in eachindex(x)
        if !(x[i]>BC.X1 || x[i]<BC.X0)
        else
            v[i] *= -1
        end
    end
end

@inline function delta(dX, BC::Union{ReflectiveBoundary,OpenBoundary})
    return dX
end

@inline function delta(dX, BC::PeriodicBoundary)
    if 2*dX > BC.L
        dX -= BC.L
    elseif 2*dX < -BC.L
        dX += BC.L
    end
    dX
end


@inline function distance(a, b, BC::GeneralBC)
    dX = delta(a[1]-b[1],BC.X)
    dY = delta(a[2]-b[2],BC.Y)
    dZ = delta(a[3]-b[3],BC.Z)
    r2 = dX*dX+dY*dY+dZ*dZ
    r = sqrt(r2)
    dr = (dX, dY, dZ)
    return r, r2, dr
end

function lattice(BC::SimulationBoundaries)
    throw("lattice fuction not defined for type $(typeof(BC)). Define it as lattice(BC::$(typeof(BC)))")
end

function lattice(BC::GeneralBC)
    f(BC) = BC.X.L, 0.0, 0.0, 0.0, BC.Y.L, 0.0, 0.0, 0.0, BC.Z.L
    str1 = "Lattice=\""
    for i in f(BC)
        str1 = "$str1$i "
    end
    str1 = "$(str1[1:end-1])\""
    return "$str1 Origin=\"$(BC.X.X0) $(BC.Y.X0) $(BC.Z.X0)\""
end



#
