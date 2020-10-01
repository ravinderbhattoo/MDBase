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

function write_trajectory_xyz(filename, res::ODESolution, params)
    mkpath("$(filename)_ovito_files/")
    for i in 1:length(res.t)
        t = res.t[i]
        data = res(t)
        v = data.x[1]'
        u = data.x[2]'
        N = params.S.N
        out = Array{Any,2}(undef, N, 11)
        out[:,6:8] = v
        out[:,3:5] = u
        out[:,9:11] = params.S.acc[:,:,i]'
        out[:,1] = params.S.sim.a_ids
        out[:,2] = params.S.sim.m_ids
        lattice_ = lattice(params.S.sim.boundary_condition)
        comment = "$lattice_ Properties=species:S:1:molecule:S:1:pos:R:3:Velocity:R:3:Acceleration:R:3"
        writedata("$(filename)_ovito_files/frame_$(string(Int(1e5+i))[2:end]).data", out, comment=comment)
    end
    print("Trajectory has been written.")
end

function write_trajectory_pvt(filename, res::ODESolution, params)
    pvd = paraview_collection(filename)
    mkpath("$(filename)_pvd_files/")
    for i in 1:length(res.t)
        t = res.t[i]
        data = res(t)
        v = data.x[1]'
        x = data.x[2]'
        N = Int(size(data, 1)/2)
        a = params.S.acc[:,:,i]'
        vtkfile = vtk_grid("$(filename)_pvd_files/frame_$i", x[:,1], x[:,2], x[:,3], MeshCell[])
        vtkfile["Velocity"] = (v[:,1], v[:,2], v[:,3])
        vtkfile["Acceleration"] = (a[:,1], a[:,2], a[:,3])
        pvd[i] = vtkfile
    end
    vtk_save(pvd)
    print("Trajectory has been written.")
end
