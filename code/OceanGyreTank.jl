using Plots
using Oceananigans
using Oceananigans.Units

p = (
    cᴰ = 2.5e-3, # dimensionless drag coefficient
    ρₐ = 1.225,  # kg m⁻³, average density of air at sea-level
    ρₒ = 1026,   # kg m⁻³, average density at the surface of the world ocean
    L = 0.23meters,
    H = 0.15meters,
    Ny = 100,
    Nx = 100,
    Nz = 20,
)


grid = RectilinearGrid(
    size=(p.Nx, p.Ny, p.Nz),
    x=(-p.L, p.L),
    y=(-p.L, p.L),
    z=(-p.H,0),
    topology=(Bounded, Bounded, Bounded)
)


radius(x,y) = sqrt(x^2 + y^2)
U(x,y,L)    = sin(π * radius(x,y) / L)*exp(1im*angle(x+y*1im))
Qᵘ(x,y,z,t,p) = radius(x,y)<p.L ?   imag(p.ρₐ / p.ρₒ * p.cᴰ * U(x,y,p.L) * abs(U(x,y,p.L))) : 0 # m² s⁻²
Qᵛ(x,y,z,t,p) = radius(x,y)<p.L ? - real(p.ρₐ / p.ρₒ * p.cᴰ * U(x,y,p.L) * abs(U(x,y,p.L))) : 0 # m² s⁻²

u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵘ, parameters=p))
v_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵛ, parameters=p))

coriolis = BetaPlane(f₀=2.0944,β=1.6)


model = NonhydrostaticModel(grid = grid,
                            tracers = (:T, :S),
                            coriolis = coriolis,
                            buoyancy = SeawaterBuoyancy(),
                            closure = ScalarDiffusivity(ν=1e-6,κ=1e-6),
                            boundary_conditions = (u=u_bcs,v=v_bcs))

set!(model,T=25,S=0)

simulation = Simulation(model, Δt = 0.1second, stop_time = 10seconds)

simulation.output_writers[:velocities] =
    JLD2OutputWriter(model, model.velocities, prefix = "../data/OceanGyreTank",
                     schedule=TimeInterval(0.1second), force = true)

run!(simulation)

u = FieldTimeSeries("../data/OceanGyreTank.jld2", "u")
v = FieldTimeSeries("../data/OceanGyreTank.jld2", "v")

xu,yu,zu = nodes(u)
xv,yv,zv = nodes(v)


anim = @animate for (i, t) in enumerate(u.times)
    p1 = contour(xu,yu,u[:,:,1,i]'.*1e3; xlabel="x [m]", ylabel="y [m]", clim = (-0.25,0.25), fill=true, linewidth=0.2)
    p2 = contour(xv,yv,v[:,:,1,i]'.*1e3; xlabel="x [m]", ylabel="y [m]", clim = (-0.25,0.25), fill=true, linewidth=0.2)
    plot(p1,p2)
end

gif(anim, "../img/animation.gif", fps = 10)