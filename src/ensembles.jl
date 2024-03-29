# exports
export Ensemble
export ENS
export ddu, ddu!, apply!

abstract type Ensemble end

#######################################################
# >>>>>>> ENS                                         #
#######################################################

struct ENS <: Ensemble
    etypes::Array{EoMType,1}
    ctypes::Array{CallBackType,1}
end

function Base.show(stream::IO, ENS::Ensemble)
    println(stream, "Ensemble:")
    println(stream, "  a) EoM types:")
    ind = 1
    for type in ENS.etypes
        println(stream, "    $ind.\t$type")
        ind += 1
    end
    println(stream, "  b) CallBack types:")
    ind = 1
    for type in ENS.ctypes
        println(stream, "    $ind.\t$type")
        ind += 1
    end
end

function ENS(ectypes::Array{<:Union{EoMType, CallBackType},1})
    ctypes = []
    etypes = []
    for arg in ectypes
        if typeof(arg) <: CallBackType
            push!(ctypes, arg)
        elseif typeof(arg) <: EoMType
            push!(etypes, arg)
        else
            throw("Not a valid parameter for ensemble. [$arg] It should be of type: $(Union{EoMType, CallBackType})")
        end
    end
    ENS(etypes,ctypes)
end


@inline function ddu(v, u, params, t, ensemble::ENS)
    ddu = 0*v
    for etype in ensemble.etypes
        ddu!(ddu, v, u, params, t, etype)
    end
    return ddu
end

@inline function ddu!(ddu, v, u, params, t, ensemble::ENS)
    for etype in ensemble.etypes
        ddu!(ddu, v, u, params, t, etype)
    end
end

@inline function apply!(v, u, params, t, ensemble::E) where E <: Ensemble
    for ctype in ensemble.ctypes
        apply!(v, u, params, t, ctype)
    end
end

#######################################################
# <<<<<<< ENS                                         #
#######################################################
