# FIGURE 2 plots S ~ eps for mu = 2,10 and xm = .1,1,2 and κm = 20,200
# (variables are named f1 because it's the first figure of results)

include("larva-base.jl")

f1_p0 = copy(p0) # base parameter set
f1_Lϵ = 50 # number of ϵ values to use
f1_pars = OrderedDict( # parameters to be varied
    "ϵ" => collect(range(0.,1,length = f1_Lϵ)),
    "χm" => [.1,1,2],
    "κm" => [20.,200.],
    "μ" => [2.,10.]
)

# do the parameter sweep, and time it
f1_time = time()
f1_df = parsweep_multi(larva_fates,f1_p0,f1_pars,[:S,:D,:W],other_args)
f1_time = time() - f1_time

# make figures
f1a = Gadfly.plot(
    filter(raw -> (raw.μ < 5),f1_df),
    y = :S,
    Scale.y_log10,
    x = :ϵ,
    linestyle = :χm,
    color = :κm,
    Scale.color_discrete_manual(my_colors[1:2]...), #(colormap = Colors.colormap("Grays",logscale = true)),
    Geom.line,
    Coord.cartesian(xmin = 0,xmax = 1,ymin = -3.5,ymax = 0),
    Guide.xlabel("ε"),
    Guide.ColorKey(title = "κ")
)

f1b = Gadfly.plot(
    filter(raw -> (raw.μ > 5),f1_df),
    y = :S,
    Scale.y_log10,
    x = :ϵ,
    linestyle = :χm,
    color = :κm,
    Scale.color_discrete_manual(my_colors[1:2]...), #(colormap = Colors.colormap("Grays",logscale = true)),
    Geom.line,
    Coord.cartesian(xmin = 0,xmax = 1,ymin = -3.5,ymax = 0),
    Guide.xlabel("ε"),
    Guide.ColorKey(title = "κ")
)
