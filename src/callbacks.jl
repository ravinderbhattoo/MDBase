export mdcallbackset, get_view_from_u, set_u

@inline function get_view_from_u(u::T) where T <: AbstractArray
    vx = @view u.x[1][1,:]
    vy = @view u.x[1][2,:]
    vz = @view u.x[1][3,:]
    ux = @view u.x[2][1,:]
    uy = @view u.x[2][2,:]
    uz = @view u.x[2][3,:]
    return ux, uy, uz, vx, vy, vz
end

@inline function get_view_from_u(u::Array{T,2}) where T <: AbstractFloat
    ux = @view u[1,:]
    uy = @view u[2,:]
    uz = @view u[3,:]
    return ux, uy, uz
end


function mdcallbackset(;list=false)
    cb1 = DiscreteCallback(updates_condition_f, updates_affect_f!, save_positions=(false,false))
    cb2 = DiscreteCallback(saveacc_condition_f, saveacc_affect_f!, save_positions=(false,false))
    if !list
        CallbackSet(cb1, cb2)
    else
        return [cb1, cb2]
    end
end

function mdcallbackset(cbs)
    lisT = mdcallbackset(list=true)
    for i in cbs
        push!(lisT, i)
    end
    CallbackSet(lisT...)
end

############################################################
# >>>>> Updates                                            #
############################################################
function updates_condition_f(u,t,integrator)
    true
end

function updates_affect_f!(integrator)
    ensemble, params = integrator.p.ensemble, integrator.p.params
    params.M.step += 1
    N = params.S.N
    ux, uy, uz, vx, vy, vz = get_view_from_u(integrator.u)

    # apply boundary conditions
    apply_simulation_bc!(ux, uy, uz, vx, vy, vz, params.S.sim.boundary_condition)

    # apply temperature updation required for some tasks
    params.M.ke = 0.5sum(@. params.S.sim.mass*(vx^2 + vy^2 + vz^2))
    params.M.Temperature = 2params.M.ke/(3params.S.N*params.S.kb)

    # updating u and v
    integrator.sol.u[end].x[1] .= integrator.u.x[1]
    integrator.sol.u[end].x[2] .= integrator.u.x[2]
end
############################################################
# <<<<< Updates                                            #
############################################################

############################################################
# >>>>> Save acceleration                                  #
############################################################
function saveacc_condition_f(u,t,integrator)
    ensemble, params = integrator.p.ensemble, integrator.p.params
    params.M.step%params.S.sim.save_every==0 || params.M.step==1
end
function saveacc_affect_f!(integrator)
    ensemble, params = integrator.p.ensemble, integrator.p.params
    du = get_du(integrator)
    params.S.acc[:,:,params.M.step] .= du.x[1]
end
############################################################
# <<<<< Save acceleration                                  #
############################################################
