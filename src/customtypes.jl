export CircularRange

struct CircularRange
    start::Int64
    stop::Int64
    length::Int64
end

function CircularRange(start, stop)
    CircularRange(start, stop, stop-start+1)
end

function CircularRange(stop)
    CircularRange(1, stop, stop)
end


Base.iterate(m::CircularRange, state=1) = state > m.length  ? nothing : (m.start + state - 1, state+1)

Base.length(m::CircularRange) = m.length
Base.lastindex(m::CircularRange) = m.stop
Base.firstindex(m::CircularRange) = m.start

Base.getindex(m::CircularRange, i::Int64) = m.start + (m.length +  (i-1)%m.length)%m.length

Base.show(stream::IO, m::CircularRange) where S = println(stream, "$(m.start):$(m.stop)")
