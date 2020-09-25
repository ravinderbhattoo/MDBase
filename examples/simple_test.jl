using MDBase

MDBase.States.FreeRun = true

u0 = 10*rand(Types.F,3,1_000)
v0 = 0.01*u0

mass = ones(Types.F, size(u0, 2))
pbc = CubicPBC(Types.F(10.0))
potentials = [Dummy()]
ensembles = [ENS([],[])]

sim = MDSim(u0, v0, mass, potentials, pbc, Δτ=0.1)

res = simulate(sim, 100, ensembles, verbose=true)

write_trajectory_pvt("./output/sample", res)
write_trajectory_xyz("./output/sample", res)



#
