macro def(name, definition)
    return quote
        macro $(esc(name))()
            esc($(Expr(:quote, definition)))
        end
    end
end

# =============================================
# export Blocks
#
# struct Block
#     x::Array{Float64, 3, N} where N
#     v::Array{Int64, 3, N} where N
#     a_ids::Array{Int64, 1, N} where N
# end
#
# struct Box
#     b::Array{Block, n1, n2, n3} where {n1, n2, n3}
#     size::Vector{Int, 3}
# end
#
#
# function neighbours(box::Box, id::vector{Int, 3})
#     a = Vector{Vector{Int, 3}}()
#     for i in id[1]-1:id[1]+1
#         for j in id[2]-1:id[2]+1
#             for k in id[3]-1:id[3]+1
#
#             end
#         end
#     end
# end
#
# =============================================



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


Base.iterate(m::CircularRange, state=1) = state > m.length  ? nothing : (m.start + state - 1, state + 1)

Base.length(m::CircularRange) = m.length
Base.lastindex(m::CircularRange) = m.stop
Base.firstindex(m::CircularRange) = m.start

Base.getindex(m::CircularRange, i::Int64) = m.start + (m.length +  (i-1)%m.length)%m.length

Base.show(stream::IO, m::CircularRange) where S = println(stream, "$(m.start):$(m.stop)")
