using Random
using Printf
using JLD2

using Oceananigans
using Oceananigans.Units: minute, minutes, hour, hours

# ## The grid
#
const H = 30
grid = RectilinearGrid(GPU();
                       size = (256, 256, 288),
                       halo = (3, 3, 3),
                          x = (0, 2π*H),
                          y = (0, 2π*H),
                          z = (-H, 0))

#
# ## Boundary conditions

# surface boundary conditions
#
u₁₀ = 8.0    # m s⁻¹, average wind velocity 10 meters above the ocean
cᴰ = 1.25e-3 # dimensionless drag coefficient
ρₐ = 1.225  # kg m⁻³, average density of air at sea-level
ρₒ = 1026.0 # kg m⁻³, average density at the surface of the world ocean
Qᵘ = - ρₐ / ρₒ * cᴰ * u₁₀ * abs(u₁₀) # m² s⁻²

# bottom drag
z₀ = 0.01 # m (roughness length)
κ = 0.4 # von Karman constant
z₁ = grid.Δzᵃᵃᶠ/2 # half the grid cell thickness
cᴰ = (κ / log(z₁ / z₀))^2 # Drag coefficient

@inline drag_u(x, y, t, u, v, p) = - p.cᴰ * √(u^2 + v^2) * u
@inline drag_v(x, y, t, u, v, p) = - p.cᴰ * √(u^2 + v^2) * v

drag_bc_u = FluxBoundaryCondition(drag_u, field_dependencies=(:u, :v), parameters=(; cᴰ))
drag_bc_v = FluxBoundaryCondition(drag_v, field_dependencies=(:u, :v), parameters=(; cᴰ))

u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵘ), bottom = drag_bc_u)
v_bcs = FieldBoundaryConditions(bottom = drag_bc_v)

# ## Coriolis
f₀ = 1e-4
coriolis = FPlane(f=f₀)

# ## Background current
#
U∞ = -0.25 # m s⁻¹

# ## Stokes drift
using Oceananigans.BuoyancyModels: g_Earth

const amplitude = 1 # m
const wavelength = 60 # m
const wavenumber = 2π / wavelength # m⁻¹
const frequency = sqrt(g_Earth * wavenumber * tanh(wavenumber * H)) # s⁻¹

## Stokes drift velocity at the surface
const Uˢ = amplitude^2 * wavenumber * frequency # m s⁻¹

# The vertical derivative of the Stokes drift is
∂z_uˢ(z, t) = Uˢ * wavenumber * sinh(2 * wavenumber * (z + H)) / (sinh(wavenumber * H))^2

# ## Forcing
#
# pressure gradient force due Coriolis in geostrophic balance
#

const pgrad = f₀ * U∞
@inline pressure_gradient(x, y, z, t) = pgrad
v_forcing = Forcing(pressure_gradient)

# ## Model instantiation
#

model = NonhydrostaticModel(; grid,
                advection = WENOFifthOrder(),
              timestepper = :RungeKutta3,
                  tracers = (:b, :Cs, :Cb),
                 coriolis = coriolis,
                 buoyancy = BuoyancyTracer(),
                  closure = AnisotropicMinimumDissipation(),
             stokes_drift = UniformStokesDrift(∂z_uˢ=∂z_uˢ),
      boundary_conditions = (; u=u_bcs, v=v_bcs),
                  forcing = (; v=v_forcing))

# ## Initial conditions
#
N² = 1.962e-4 # s⁻²
b0(z) = N² * (H+z)

## Random noise damped at top and bottom
Ξ(z) = randn() * z / model.grid.Lz * (1 + z / model.grid.Lz) # noise

## Buoyancy initial condition: a stable density gradient with random noise superposed
bᵢ(x, y, z) = b0(z) + N² * model.grid.Lz * 1e-6 * Ξ(z)

## Horizontal velocity initial condition: geostrophic flow with random noise superposed
uᵢ(x, y, z) = U∞ + 1e-5 * Ξ(z)

## Vertical velocity initial condition: random noise scaled by the friction velocity.
wᵢ(x, y, z) = sqrt(abs(Qᵘ)) * 1e-3 * Ξ(z)

## Tracer initial condition: surface tracer and bottom tracer
h = 0.5 # m
Csᵢ(x, y, z) = z < -h ? 0.0 : 1.0
Cbᵢ(x, y, z) = z > -H+h ? 0.0 : 1.0

## `set!` the `model` fields using functions or constants:
set!(model, u=uᵢ, w=wᵢ, b=bᵢ, Cs=Csᵢ, Cb=Cbᵢ)

# ## Setting up a simulation
#

simulation = Simulation(model, Δt=0.763615, stop_time=384hours)

# The `TimeStepWizard` helps ensure stable time-stepping
# with a Courant-Freidrichs-Lewy (CFL) number of 1.0.

wizard = TimeStepWizard(cfl=1.0, max_change=1.1, max_Δt=1minute)

simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

# Nice progress messaging is helpful:

## Print a progress message
progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, max(|w|) = %.1e ms⁻¹, wall time: %s\n",
                                iteration(sim),
                                prettytime(sim),
                                prettytime(sim.Δt),
                                maximum(abs, sim.model.velocities.w),
                                prettytime(sim.run_wall_time))

simulation.callbacks[:progress] = Callback(progress_message, IterationInterval(100))

# We then set up the simulation:

# ## Output
#
## Create a NamedTuple with eddy viscosity
eddy_viscosity = (; νₑ = model.diffusivity_fields.νₑ)

simulation.output_writers[:slices_xy] =
    JLD2OutputWriter(model, merge(model.velocities, model.tracers, eddy_viscosity),
                         filename = "slices_xy.jld2",
                         indices = (:, :, grid.Nz-3),
                         schedule = TimeInterval(3minute),
               overwrite_existing = false)

simulation.output_writers[:slices_xy2] =
    JLD2OutputWriter(model, merge(model.velocities, model.tracers, eddy_viscosity),
                         filename = "slices_xy2.jld2",
                         indices = (:, :, 3),
                         schedule = TimeInterval(3minute),
               overwrite_existing = false)

simulation.output_writers[:slices_xy3] =
    JLD2OutputWriter(model, merge(model.velocities, model.tracers, eddy_viscosity),
                         filename = "slices_xy3.jld2",
                         indices = (:, :, 48),
                         schedule = TimeInterval(3minute),
               overwrite_existing = false)

simulation.output_writers[:slices_xy4] =
    JLD2OutputWriter(model, merge(model.velocities, model.tracers, eddy_viscosity),
                         filename = "slices_xy4.jld2",
                         indices = (:, :, 53),
                         schedule = TimeInterval(3minute),
               overwrite_existing = false)

simulation.output_writers[:slices_xy5] =
    JLD2OutputWriter(model, merge(model.velocities, model.tracers, eddy_viscosity),
                         filename = "slices_xy5.jld2",
                         indices = (:, :, 58),
                         schedule = TimeInterval(3minute),
               overwrite_existing = false)

simulation.output_writers[:slices_xz] =
    JLD2OutputWriter(model, merge(model.velocities, model.tracers, eddy_viscosity),
                         filename = "slices_xz.jld2",
                         indices = (:, 1, :),
                         schedule = TimeInterval(3minute),
               overwrite_existing = false)

simulation.output_writers[:slices_yz] =
    JLD2OutputWriter(model, merge(model.velocities, model.tracers, eddy_viscosity),
                         filename = "slices_yz.jld2",
                         indices = (1, :, :),
                         schedule = TimeInterval(3minute),
               overwrite_existing = false)

fields_to_output = merge(model.velocities, model.tracers, (νₑ=model.diffusivity_fields.νₑ,))
simulation.output_writers[:fields] =
    JLD2OutputWriter(model, merge(model.velocities, model.tracers, eddy_viscosity),
                     filename = "fields.jld2",
                     schedule = TimeInterval(3hours),
           overwrite_existing = false)

u, v, w = model.velocities
U = Average(u, dims=(1, 2))
V = Average(v, dims=(1, 2))
b = Average(model.tracers.b, dims=(1, 2))
wb = Average(w * model.tracers.b, dims=(1, 2))
wu = Average(w * u, dims=(1, 2))
wv = Average(w * v, dims=(1, 2))
ww = Average(w * w, dims=(1, 2))
uu = Average((u-Field(U)) * (u-Field(U)), dims=(1, 2))
vv = Average((v-Field(V)) * (v-Field(V)), dims=(1, 2))
uv = Average((u-Field(U)) * (v-Field(V)), dims=(1, 2))
w3 = Average(w^3, dims=(1, 2))
bb = Average((model.tracers.b-Field(b)) * (model.tracers.b-Field(b)), dims=(1, 2))
wbsb = Average(-∂z(model.tracers.b) * model.diffusivity_fields.νₑ, dims=(1, 2))
wusb = Average(-∂z(u) * model.diffusivity_fields.νₑ, dims=(1, 2))
wvsb = Average(-∂z(v) * model.diffusivity_fields.νₑ, dims=(1, 2))
Cs = Average(model.tracers.Cs, dims=(1, 2))
wCs = Average(w * model.tracers.Cs, dims=(1, 2))
wCssb = Average(-∂z(model.tracers.Cs) * model.diffusivity_fields.νₑ, dims=(1, 2))
Cb = Average(model.tracers.Cb, dims=(1, 2))
wCb = Average(w * model.tracers.Cb, dims=(1, 2))
wCbsb = Average(-∂z(model.tracers.Cb) * model.diffusivity_fields.νₑ, dims=(1, 2))

simulation.output_writers[:averages] =
    JLD2OutputWriter(model, (u=U, v=V, b=b, wb=wb, wu=wu, wv=wv, ww=ww, uu=uu, vv=vv, uv=uv, w3=w3, bb=bb, wbsb=wbsb, wusb=wusb, wvsb=wvsb,
                             Cs=Cs, wCs=wCs, wCssb=wCssb, Cb=Cb, wCb=wCb, wCbsb=wCbsb),
                     schedule = TimeInterval(3minutes),
                     filename = "averages.jld2",
           overwrite_existing = false)

checkpointer = Checkpointer(model,
                            schedule = TimeInterval(48hours),
                  overwrite_existing = false)
simulation.output_writers[:checkpointer] = checkpointer

# ## Run the model
run!(simulation, pickup=true)

