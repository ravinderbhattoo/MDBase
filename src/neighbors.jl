# exports
export find_neighbors, find_neighbors!, distance


"""
    neighbor_cells(i, j, k, N)

Calculate all neighboring cells for i, j, k cell.
"""
function neighbor_cells(i, j, k, N)
    a = Vector{Int}()
    Ni = CircularRange(N[1])
    Nj = CircularRange(N[2])
    Nk = CircularRange(N[3])
    for kk in k-1:k+1
        for jj in j-1:j+1
            for ii in i-1:i+1
                push!(a, cell_number(Ni[ii], Nj[jj], Nk[kk], N))
            end
        end
    end
    return unique(a)
end

"""
    cell_number(i, j, k, N)

Calculates cell number for given i, j, k cell indices (when in a list).
"""
function cell_number(i, j, k, N)
    return i+(j-1)*N[1]+(k-1)*N[1]*N[2]
end

"""
    get_cells(x::Array{Float64, 2}, horizon::Float64)

Fill cells with particles.
"""
function get_cells(x::Array{Float64, 2}, horizon::Float64)
    _min = minimum(x, dims=2)
    _max = maximum(x, dims=2)
    N = max.(Int.(1 .+ floor.((_max-_min)/horizon.-1.0e-6)),1)
    cells = [Vector{Int}() for i in 1:prod(N)]
    cell_neighs = Vector{Vector{Int}}(undef, prod(N))
    for k in 1:N[3]
        for j in 1:N[2]
            for i in 1:N[1]
                cell_neighs[cell_number(i, j, k, N)] = neighbor_cells(i, j, k, N)
            end
        end
    end
    for i in 1:size(x, 2)
        ii, jj, kk = max.(1, Int.(1 .+ floor.((x[:, i].-_min)/horizon.-1.0e-6)))
        push!(cells[cell_number(ii, jj, kk, N)], i)
    end
    return cells, cell_neighs
end


"""
    find_neighbors!(neighbors::Array{Int64, 2}, x::Array{Float64, 2}, horizon::Float64)

Calculate neighbors for each particle (inplace).
"""
function find_neighbors!(neighbors::Array{Int64, 2}, x::Array{Float64, 2}, horizon::Float64, BC)
    max_neighs = size(neighbors, 1)
    h2 = horizon^2
    cells, cell_neighs = get_cells(x, horizon)
    for cell_i in 1:length(cells)
        for ca in cells[cell_i]
            ind = 1
            for neighs in cell_neighs[cell_i]
                for fa in cells[neighs]
                    if ca!=fa
                        dX = delta(x[1, fa] - x[1, ca], BC.X)
                        dY = delta(x[2, fa] - x[2, ca], BC.Y)
                        dZ = delta(x[3, fa] - x[3, ca], BC.Z)
                        r2 = dX*dX+dY*dY+dZ*dZ
                        if r2<h2
                            ind += 1
                            if ind<=max_neighs
                                neighbors[ind, ca] = fa
                            else
                                error("Neighbors are more($ind) than capacity($max_neighs) allocated. Please increase neighbors maximum capacity.")
                            end
                        end
                    end
                end
            end
        end
    end
    neighbors[1,:] = sum(neighbors.!=0,dims=1)
    nothing
end


"""
    find_neighbors(x::Array{Float64, 2}, horizon::Float64)

Calculate neighbors for each particle.
"""
function find_neighbors(x::Array{Float64, 2}, horizon::Float64, BC, max_neigh::Int64; hard_max::Int64=0)::Array{Int64, 2}
    max_neigh = max(hard_max, min(max_neigh, size(x, 2)))
    neighbors = zeros(Int64, max_neigh, size(x, 2))
    find_neighbors!(neighbors, x, horizon, BC)
    return neighbors
end
