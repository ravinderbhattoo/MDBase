using MDBase

MDBase.States.FreeRun = true

u0 = 10*rand(Types.F,3,1_000)
v0 = 0.5 .- rand(Types.F,3,1_000)

mass = ones(Types.F, size(u0, 2))
pbc = CubicPBC(Types.F(10.0))
potentials = [Dummy()]
ensembles = [ENS([],[])]

sim = MDSim(u0, v0, mass, potentials, pbc, Δτ=0.1, save_every=1, thermo_save_every=1)

res, parameters = simulate(100, sim, ensembles, verbose=true)

# prob, dt, saveat, parameters = _stage(sim, 100, ensembles, verbose=true)
# using BenchmarkTools
# @time MDBase._simulate(prob, dt, saveat, parameters)

write_trajectory_pvt("./output/sample", res, parameters)
write_trajectory_xyz("./output/sample", res, parameters)


#
