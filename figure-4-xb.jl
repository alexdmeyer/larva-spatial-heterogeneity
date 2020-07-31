## FIGURE 4: a-c are S,W,D ~ xb with eps = 1, different mu's
# d-i are S,W,D ~ xb with eps = .1, different xm's, mu = 2 or 10

include("larva-base.jl")

f3u_p0 = copy(p0) # base parameters for UNIFORM mortality
f3h_p0 = copy(p0); f3h_p0["ϵ"] = .1; # base for HETEROGENEOUS mortality
f3_Lxb = 50 # number of CBL widths χb to use

# low mortality cases require larger spatial domain since larvae survive
# to travel further out. We run those separately because large domain
# unnecessary for higher mortality cases
f3u_pars_lowMortality = OrderedDict(
    "χb" => vcat(0.,range(.5,10,length = f3_Lxb)), #collect(range(0.1,10.,length = f3a_Lxb)),
    "L" => [100.],
    "μ" => [0.,2.]
)
f3u_pars = OrderedDict(
    "χb" => vcat(0.,range(.5,10,length = f3_Lxb)), #collect(range(0.1,10.,length = f3a_Lxb)),
    "μ" => [10.,],
    "L" => [50.],
)
f3h_pars_lowMortality = OrderedDict(
    "χb" => vcat(0.,range(.5,10,length = f3_Lxb)), #collect(range(0.1,10.,length = f3a_Lxb)),
    "χm" => [.1,1,2],
    "L" => [100.],
    "μ" => [2.]
)
f3h_pars = OrderedDict(
    "χb" => vcat(0.,range(.5,10,length = f3_Lxb)), #collect(range(0.1,10.,length = f3a_Lxb)),
    "L" => [50.],
    "χm" => [.1,1.,2.],
    "μ" => [10.]
)

# do the parameter sweeps, and time them all together
f3_time = time()
f3u_df = parsweep_multi(larva_fates,f3u_p0,f3u_pars,[:S,:D,:W],other_args)
f3u_df2 = parsweep_multi(larva_fates,f3u_p0,f3u_pars_lowMortality,[:S,:D,:W],other_args)
f3h_df = parsweep_multi(larva_fates,f3h_p0,f3h_pars,[:S,:D,:W],other_args)
f3h_df2 = parsweep_multi(larva_fates,f3h_p0,f3h_pars_lowMortality,[:S,:D,:W],other_args)
f3_time = time() - f3_time

# make figures
f3a = Gadfly.plot(
    filter(row -> (row.μ < ρ),f3u_df2),
    layer(
        x = :χb,
        y = :S,
        Theme(default_color = my_colors[1]),
        Geom.line
        ),
    layer(
        x = :χb,
        y = :W,
        Theme(default_color = my_colors[2]),
        Geom.line
        ),
    Coord.cartesian(ymin = -2.5,ymax = 0),
    Scale.y_log10,
    Guide.xlabel("xb"),
    Guide.ylabel("probability"),
    Guide.ManualColorKey("fate",["S","W"],my_colors[1:3])
)

f3b = Gadfly.plot(
    filter(row -> abs(row.μ-2) < ρ,f3u_df2),
    layer(
        x = :χb,
        y = :S,
        Theme(default_color = my_colors[1]),
        Geom.line
        ),
    layer(
        x = :χb,
        y = :W,
        Theme(default_color = my_colors[2]),
        Geom.line
        ),
    layer(
        x = :χb,
        y = :D,
        Theme(default_color = my_colors[3]),
        Geom.line
        ),
    Scale.y_log10,
    Coord.cartesian(ymin = -2.5,ymax = 0),
    Guide.xlabel("xb"),
    Guide.ylabel("probability"),
    Guide.ManualColorKey("fate",["S","W","M"],my_colors[1:3])
)

f3c = Gadfly.plot(
    filter(row -> abs(row.μ-10) < ρ,f3u_df),
    layer(
        x = :χb,
        y = :S,
        Theme(default_color = my_colors[1]),
        Geom.line
        ),
    layer(
        x = :χb,
        y = :W,
        Theme(default_color = my_colors[2]),
        Geom.line
        ),
    layer(
        x = :χb,
        y = :D,
        Theme(default_color = my_colors[3]),
        Geom.line
        ),
    Scale.y_log10,
    Coord.cartesian(ymin = -4,ymax = 0),
    Guide.xlabel("xb"),
    Guide.ylabel("probability"),
    Guide.ManualColorKey("fate",["S","W","M"],my_colors[1:3])
)

f3d = Gadfly.plot(
    filter(row -> (abs(row.χm-.1) < ρ),f3h_df2),
    layer(
        x = :χb,
        y = :S,
        Theme(default_color = my_colors[1]),
        Geom.line
        ),
    layer(
        x = :χb,
        y = :W,
        Theme(default_color = my_colors[2]),
        Geom.line
        ),
    layer(
        x = :χb,
        y = :D,
        Theme(default_color = my_colors[3]),
        Geom.line
        ),
    Coord.cartesian(ymin = -2.5,ymax = 0),
    Scale.y_log10,
    Guide.xlabel("xb"),
    Guide.ylabel("probability"),
    Guide.ManualColorKey("fate",["S","W","M"],my_colors[1:3])
)

f3e = Gadfly.plot(
    filter(row -> (abs(row.χm - 1) < ρ),f3h_df2),
    layer(
        x = :χb,
        y = :S,
        Theme(default_color = my_colors[1]),
        Geom.line
        ),
    layer(
        x = :χb,
        y = :W,
        Theme(default_color = my_colors[2]),
        Geom.line
        ),
    layer(
        x = :χb,
        y = :D,
        Theme(default_color = my_colors[3]),
        Geom.line
        ),
    Scale.y_log10,
    Coord.cartesian(ymin = -2.5,ymax = 0),
    Guide.xlabel("xb"),
    Guide.ylabel("probability"),
    Guide.ManualColorKey("fate",["S","W","M"],my_colors[1:3])
)

f3f = Gadfly.plot(
    filter(row -> (abs(row.μ-2) < ρ) & (abs(row.χm - 2) < ρ),f3h_df2),
    layer(
        x = :χb,
        y = :S,
        Theme(default_color = my_colors[1]),
        Geom.line
        ),
    layer(
        x = :χb,
        y = :W,
        Theme(default_color = my_colors[2]),
        Geom.line
        ),
    layer(
        x = :χb,
        y = :D,
        Theme(default_color = my_colors[3]),
        Geom.line
        ),
    Coord.cartesian(ymin = -2.5,ymax = 0),
    Scale.y_log10,
    Guide.xlabel("xb"),
    Guide.ylabel("probability"),
    Guide.ManualColorKey("fate",["S","W","M"],my_colors[1:3])
)

f3g = Gadfly.plot(
    filter(row -> (abs(row.χm-.1) < ρ),f3h_df),
    layer(
        x = :χb,
        y = :S,
        Theme(default_color = my_colors[1]),
        Geom.line
        ),
    layer(
        x = :χb,
        y = :W,
        Theme(default_color = my_colors[2]),
        Geom.line
        ),
    layer(
        x = :χb,
        y = :D,
        Theme(default_color = my_colors[3]),
        Geom.line
        ),
    Coord.cartesian(ymin = -2.5,ymax = 0),
    Scale.y_log10,
    Guide.xlabel("xb"),
    Guide.ylabel("probability"),
    Guide.ManualColorKey("fate",["S","W","M"],my_colors[1:3])
)

f3h = Gadfly.plot(
    filter(row -> (abs(row.χm - 1) < ρ),f3h_df),
    layer(
        x = :χb,
        y = :S,
        Theme(default_color = my_colors[1]),
        Geom.line
        ),
    layer(
        x = :χb,
        y = :W,
        Theme(default_color = my_colors[2]),
        Geom.line
        ),
    layer(
        x = :χb,
        y = :D,
        Theme(default_color = my_colors[3]),
        Geom.line
        ),
        Scale.y_log10,
    Coord.cartesian(ymin = -2.5,ymax = 0),
    Guide.xlabel("xb"),
    Guide.ylabel("probability"),
    Guide.ManualColorKey("fate",["S","W","M"],my_colors[1:3])
)

f3i = Gadfly.plot(
    filter(row -> (abs(row.χm - 2) < ρ),f3h_df),
    layer(
        x = :χb,
        y = :S,
        Theme(default_color = my_colors[1]),
        Geom.line
        ),
    layer(
        x = :χb,
        y = :W,
        Theme(default_color = my_colors[2]),
        Geom.line
        ),
    layer(
        x = :χb,
        y = :D,
        Theme(default_color = my_colors[3]),
        Geom.line,
        ),
    Scale.y_log10,
    Coord.cartesian(ymin = -2.5,ymax = 0),
    Guide.xlabel("xb"),
    Guide.ylabel("probability"),
    Guide.ManualColorKey("fate",["S","W"],my_colors[1:2])
)
