# FIGURE 3 plot S ~ κ̄ for ϵ = 0,.25,.5,1 and χm = .1,1,2 and μ = 2,10

include("larva-base.jl")

f2_p0 = copy(p0) # base parameters
f2_Lκ = 50 # number of κ values to use
f2_pars = OrderedDict(
    "ϵ" => collect([0.,.25,.5,1.]),
    "κ̄" => 10. .^range(0.,3.,length = f2_Lκ),
    "xm" => [.1,1.,2.],
    "μ" => [2.,10.]
)

# do the parameter sweep, and time it
f2_time = time()
f2_df2 = parsweep_multi(larva_fates,f2_p0,f2_pars,[:S,:M,:W],other_args)
f2_time = time() - f2_time

# make figures
f2a = Gadfly.plot(
    filter(row -> (row.μ < 5.) & (row.xm < .5),f2_df2),
    y = :S,
    x = :κ̄,
    color = :ϵ,
    Scale.color_discrete,
    Scale.color_discrete_manual(my_colors...),
    Scale.y_log10,
    Scale.x_log10,
    Geom.line,
    Coord.cartesian(ymin = -4,ymax = 0),
    Guide.xlabel("κ̄"),
    Guide.ColorKey(title = "ε"),
)

f2b = Gadfly.plot(
    filter(row -> (row.μ < 5.) & (1.5 > row.xm > .5),f2_df2),
    y = :S,
    x = :κ̄,
    color = :ϵ,
    Geom.line,
    Coord.cartesian(ymin = -4,ymax = 0),
    Scale.color_discrete_manual(my_colors...),
    Scale.y_log10,
    Scale.x_log10,
    Guide.xlabel("κ̄"),
    Guide.ColorKey(title = "ε"),
)

f2c = Gadfly.plot(
    filter(row -> (row.μ < 5.) & (row.xm > 1.5),f2_df2),
    y = :S,
    x = :κ̄,
    color = :ϵ,
    Geom.line,
    Coord.cartesian(ymin = -4,ymax = 0),
    Scale.color_discrete_manual(my_colors...),
    Scale.y_log10,
    Scale.x_log10,
    Guide.xlabel("κ̄"),
    Guide.ColorKey(title = "ε"),
)

f2d = Gadfly.plot(
    filter(row -> (row.μ > 5.) & (row.xm < .5),f2_df2),
    y = :S,
    x = :κ̄,
    color = :ϵ,
    Geom.line,
    Coord.cartesian(ymin = -4,ymax = 0),
    Scale.color_discrete_manual(my_colors...),
    Scale.y_log10,
    Scale.x_log10,
    Guide.xlabel("κ̄"),
    Guide.ColorKey(title = "ε"),
)

f2e = Gadfly.plot(
    filter(row -> (row.μ > 5.) & (1.5 > row.xm > .5),f2_df2),
    y = :S,
    x = :κ̄,
    color = :ϵ,
    Geom.line,
    Coord.cartesian(ymin = -4,ymax = 0),
    Scale.color_discrete_manual(my_colors...),
    Scale.y_log10,
    Scale.x_log10,
    Guide.xlabel("κ̄"),
    Guide.ColorKey(title = "ε"),
)

f2f = Gadfly.plot(
    filter(row -> (row.μ > 5.) & (row.xm > 1.5),f2_df2),
    y = :S,
    x = :κ̄,
    color = :ϵ,
    Geom.line,
    Coord.cartesian(ymin = -4,ymax = 0),
    Scale.color_discrete_manual(my_colors...),
    Scale.y_log10,
    Scale.x_log10,
    Guide.xlabel("κ̄"),
    Guide.ColorKey(title = "ε"),
)
