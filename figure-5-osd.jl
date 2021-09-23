# Figure 5: density plots of offshore duration θ for μ = 2,10 and
# ϵ = 0.1, 1.0. This code takes a very long time to run

include("larva-base.jl")
p0["κm"] = 100.

M = Int(1e4) # number of simulations
dt = .025 # time step
tmax = 1. # final time
tvec = collect(0:dt:tmax) # time mesh
DD = [1,10]

# w = weak, s = strong, u = uniform, h = heterogeneous,
# cbl = coastal boundary layer (in the end, CBL made little difference)
# cases = ["wu","wu_cbl","su","su_cbl","wh","wh_cbl","sh","sh_cbl"]
cases = ["wu","wh"]
# create a dictionary of parameter dictionaries so we can index parameter sets
# by case name
pars = Dict(case => copy(p0) for case in cases)
pars["wu"]["μ"] = 6.; pars["wu"]["ϵ"] = 1.;
pars["wh"]["μ"] = 6.; pars["wh"]["ϵ"] = 0.;

# pars["su"]["μ"] = 10.; pars["su"]["ϵ"] = 1.
# pars["sh"]["μ"] = 10.; pars["sh"]["ϵ"] = 0.

# pars["wu_cbl"]["μ"] = 2.; pars["wu_cbl"]["ϵ"] = 1.; pars["wu_cbl"]["χb"] = 5.
# pars["wh_cbl"]["μ"] = 2.; pars["wh_cbl"]["ϵ"] = 0.; pars["wh_cbl"]["χb"] = 5.
# pars["su_cbl"]["μ"] = 10.; pars["su_cbl"]["ϵ"] = 1.; pars["su_cbl"]["χb"] = 5.
# pars["sh_cbl"]["μ"] = 10.; pars["sh_cbl"]["ϵ"] = 0.; pars["sh_cbl"]["χb"] = 5.

# preallocate ordered dictionaries to store results (ordering is so indices match)
XX = OrderedDict(); # store simulated time series
df = OrderedDict(); # store dataframes of time series
Time = OrderedDict(); # time at which killing occurs

# simulate M trials, all ending in settling. Can't run sims in parallel, but
# can at least run cases in parallel
Threads.@threads for case in cases
    XX[case],df[case],Time[case] = KilledDiffusion1d_DistTraveled(M,pars[case],f,g,rand(M),tvec,BCs,k,kill_names,time_beyond_D = DD,kill_type = "settling",step_method = em1d_step)
    println("Case " * case * " complete")
    println(case * " run-time = " *string(Time[case]))
    println(case * " thread = " * string(Threads.threadid()))
end

# store everything in big dataframe
bigdf = DataFrame()
# bigdf[!,:μ] = repeat([2.,10.],inner = 2*M,outer = 2)
# bigdf[!,:ϵ] = repeat([1.,.1],inner = 4*M,)
bigdf[!,:ϵ] = repeat([1.,0],inner = M,)
# bigdf[!,:χb] = repeat([0.,5.],inner = M,outer = 4)
for D in DD
    bigdf[!,Symbol("θ"*string(D))] = vcat([df[case][:,Symbol("TimeBeyond"*string(D))] for case in cases]...)
end

bigdf[!,:Xmax] = vcat([df[case][:,:MaxDistance] for case in cases]...)

a = plot(
    layer(filter(r -> (r.ϵ < .5),bigdf),x = :θ1,Geom.histogram,Theme(default_color = my_colors[1],alphas = [.5])),
    layer(filter(r -> (r.ϵ > .5),bigdf),x = :θ1,Geom.histogram,Theme(default_color = my_colors[2],alphas = [.5])),
    Coord.cartesian(xmin = 0,xmax = 1),
    Guide.xlabel("time beyond x = 1, θ(1)"),
    Guide.ylabel("number of trials")
)
# draw(PDF(figz * "/fig-theta1.pdf",12cm,10cm),a)

b = plot(
    layer(filter(r -> (r.ϵ < .5),bigdf),x = :θ10,Geom.histogram,Theme(default_color = my_colors[2],alphas = [.5],discrete_highlight_color=c->colorant"black")),
    layer(filter(r -> (r.ϵ > .5),bigdf),x = :θ10,Geom.histogram,Theme(default_color = my_colors[1],alphas = [.5],discrete_highlight_color=c->colorant"black")),
    Coord.cartesian(xmin = 0,xmax = 1),
    Guide.xlabel("time beyond x = 10, θ(10)"),
    Guide.ylabel("settled larvae"),
    Guide.manual_color_key("ε",["1.0","0.0"],my_colors[1:2],shape = [Shape.square],size = [1.5mm]),
    latex_fonts
)

draw(PDF(figz * "/fig-theta10-final.pdf",15cm,10cm),b)

c1 = plot(filter(r -> (r.ϵ > .5),bigdf),x = :θ10,color = :ϵ, Scale.color_discrete_manual(my_colors[1]),Geom.histogram,Guide.xlabel("Time beyond x = 10"),Guide.ylabel("Settled larvae"),latex_fonts,Coord.Cartesian(ymin = 0,ymax = 1000))
c2 = plot(filter(r -> (r.ϵ < .5),bigdf),x = :θ10,color = :ϵ, Scale.color_discrete_manual(my_colors[2]),Geom.histogram,Guide.xlabel("Time beyond x = 10"),Guide.ylabel("Settled larvae"),latex_fonts,Coord.Cartesian(ymin = 0,ymax = 1000))
draw(PDF(figz*"/fig-theta10-unif.pdf",15cm,10cm),c1)
draw(PDF(figz*"/fig-theta10-het.pdf",15cm,10cm),c2)

# draw(PNG("/Users/alexandermeyer/Dropbox/2020-2021/Postdoc/UMD-Bruns/fig-hist-unif.png",15cm,10cm),c1)
# draw(PNG("/Users/alexandermeyer/Dropbox/2020-2021/Postdoc/UMD-Bruns/fig-hist-het.png",15cm,10cm),c2)
