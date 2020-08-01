# Figure 5: density plots of offshore duration θ for μ = 2,10 and
# ϵ = 0.1, 1.0. This code takes a very long time to run

include("larva-base.jl")

M = 1e4 # number of simulations
dt = .01 # time step
tmax = 1. # final time
tvec = collect(0:dt:tmax) # time mesh

# w = weak, s = strong, u = uniform, h = heterogeneous,
# cbl = coastal boundary layer (in the end, CBL made little difference)
cases = ["wu","wu_cbl","su","su_cbl","wh","wh_cbl","sh","sh_cbl"]
# create a dictionary of parameter dictionaries so we can index parameter sets
# by case name
pars = Dict(case => copy(p0) for case in cases)
pars["wu"]["μ"] = 2.; pars["wu"]["ϵ"] = 1.
pars["wu_cbl"]["μ"] = 2.; pars["wu_cbl"]["ϵ"] = 1.; pars["wu_cbl"]["xb"] = 5.
pars["su"]["μ"] = 10.; pars["su"]["ϵ"] = 1.
pars["su_cbl"]["μ"] = 10.; pars["su_cbl"]["ϵ"] = 1.; pars["su_cbl"]["xb"] = 5.
pars["wh"]["μ"] = 2.; pars["wh"]["ϵ"] = .1
pars["wh_cbl"]["μ"] = 2.; pars["wh_cbl"]["ϵ"] = .1; pars["wh_cbl"]["xb"] = 5.
pars["sh"]["μ"] = 10.; pars["sh"]["ϵ"] = .1
pars["sh_cbl"]["μ"] = 10.; pars["sh_cbl"]["ϵ"] = .1; pars["sh_cbl"]["xb"] = 5.

# preallocate ordered dictionaries to store results (ordering is so indices match)
XX = OrderedDict() # store simulated time series
df = OrderedDict() # store dataframes of time series
Time = OrderedDict() # time at which killing occurs

# simulate M trials, all ending in settling. Can't run sims in parallel, but
# can at least run cases in parallel
Threads.@threads for case in cases
    XX[case],df[case],Time[case] = KilledDiffusion1d_DistTraveled(M,pars[case],f,g,rand(M),tvec,BCs,k,kill_names,kill_type = "settling",step_method = srk1d_step)
    println("Case " * case * " complete")
    println(case * " run-time = " *string(Time[case]))
    println(case * " thread = " * string(Threads.threadid()))
end

# store everything in big dataframe
bigdf = DataFrame()
bigdf[!,:μ] = repeat([2.,10.],inner = 2*M,outer = 2)
bigdf[!,:ϵ] = repeat([1.,.1],inner = 4*M,)
bigdf[!,:xb] = repeat([0.,5.],inner = M,outer = 4)
bigdf[!,:θ] = vcat([df[case][:,:TimeBeyondHab] for case in cases]...)

# now plot the results

f4a = Gadfly.plot(
    filter(row -> (row.μ < 5) & (row.xb < 1), bigdf),
    x = :θ,
    color = :ϵ,
    Geom.density(bandwidth = .005),
    Coord.Cartesian(xmin = 0,xmax = 1,ymin = 0,ymax = 10),
    Scale.color_discrete_manual(alex_colors[1:2]...),
    Guide.xlabel("θ_OSD"),
    Guide.ColorKey(title = "ε"),
    Guide.ylabel("prob density",orientation = :vertical)
)

f4b = Gadfly.plot(
    filter(row -> (row.μ > 5) & (row.xb < 1), bigdf),
    x = :θ,
    color = :ϵ,
    Geom.density(bandwidth = .005),
    Coord.Cartesian(xmin = 0,xmax = 1,ymin = 0,ymax = 10),
    Scale.color_discrete_manual(alex_colors[1:2]...),
    Guide.xlabel("θ_OSD"),
    Guide.ColorKey(title = "ε"),
    Guide.ylabel("prob density",orientation = :vertical)
)

f4c = Gadfly.plot(
    filter(row -> (row.μ < 5) & (row.xb > 1), bigdf),
    x = :θ,
    color = :ϵ,
    Geom.density(bandwidth = .005),
    Coord.Cartesian(xmin = 0,xmax = 1,ymin = 0,ymax = 10),
    Scale.color_discrete_manual(alex_colors[1:2]...),
    Guide.xlabel("θ_OSD"),
    Guide.ColorKey(title = "ε"),
    Guide.ylabel("prob density",orientation = :vertical)
)

f4d = Gadfly.plot(
    filter(row -> (row.μ > 5) & (row.xb > 1), bigdf),
    x = :θ,
    color = :ϵ,
    Geom.density(bandwidth = .005),
    Coord.Cartesian(xmin = 0,xmax = 1,ymin = 0,ymax = 10),
    Scale.color_discrete_manual(alex_colors[1:2]...),
    Guide.xlabel("θ_OSD"),
    Guide.ColorKey(title = "ε"),
    Guide.ylabel("prob density",orientation = :vertical)
)
