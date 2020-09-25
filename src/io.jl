export write_trajectory_xyz, write_trajectory_pvt

using WriteVTK: paraview_collection, vtk_grid, MeshCell, vtk_save
using DelimitedFiles: readdlm, writedlm

function writedata(filename, vals; comment=nothing)
    N = size(vals, 1)
    open(filename; write=true) do f
        write(f, "$N\n")
        if comment!=nothing
            write(f, "$comment\n")
        else
            write(f, "\n")
        end
        writedlm(f, vals)
    end
end

function write_trajectory_xyz(filename, res::ODESolution)
    mkpath("$(filename)_ovito_files/")
    for i in 1:length(res.t)
        t = res.t[i]
        data = Array(reshape(Array(res(t)),(3,:)))'
        N = res.prob.p[2].S.N
        out = Array{Any,2}(undef, N, 11)
        out[:,6:8] = data[1:N,:]
        out[:,3:5] = data[N+1:end,:]
        out[:,9:11] = Array(reshape(res.prob.p[2].S.acc[:,i],(3,:)))'
        out[:,1] = res.prob.p[2].S.sim.a_ids
        out[:,2] = res.prob.p[2].S.sim.m_ids
        lattice_ = lattice(res.prob.p[2].S.sim.boundary_condition)
        comment = "$lattice_ Properties=species:S:1:molecule:S:1:pos:R:3:Velocity:R:3:Acceleration:R:3"
        writedata("$(filename)_ovito_files/frame_$(string(Int(1e5+i))[2:end]).data", out, comment=comment)
    end
    print("Trajectory has been written.")
end

function write_trajectory_pvt(filename, res::ODESolution)
    pvd = paraview_collection(filename)
    mkpath("$(filename)_pvd_files/")
    for i in 1:length(res.t)
        t = res.t[i]
        data = Array(reshape(Array(res(t)),(3,:)))'
        N = Int(size(data, 1)/2)
        v = data[1:N,:]
        x = data[N+1:end,:]
        a = Array(reshape(res.prob.p[2].S.acc[:,i],(3,:)))'
        vtkfile = vtk_grid("$(filename)_pvd_files/frame_$i", x[:,1], x[:,2], x[:,3], MeshCell[])
        vtkfile["Velocity"] = (v[:,1], v[:,2], v[:,3])
        vtkfile["Acceleration"] = (a[:,1], a[:,2], a[:,3])
        pvd[i] = vtkfile
    end
    vtk_save(pvd)
    print("Trajectory has been written.")
end
