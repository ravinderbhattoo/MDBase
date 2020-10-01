using Flux

struct abc
    a
    b
    c
end


a = abc(1,2,3)
b = abc(7,8,9)

p = [4,5,6]

function Base.getproperty(x::abc, f::Symbol)
    # if f==:a
    #     l = 1*p[1]
    #     return l
    # else
        getfield(x, f)
    # end
end


xs = collect(0:0.01:10)
ys = 0.1*xs

function loss()
    s = b.a*sum(@. (a.a* p[1]*xs - b.a*ys)^2)
end

loss()

ps = Flux.params(a)
gs = Flux.gradient(loss, ps)

k = collect(keys(gs.grads))

for i in k
    println(typeof(i))
end

gs.grads[k[end]]

gs.grads[k[end-1]]

#
