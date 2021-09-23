# FIGURE 2 plots S ~ eps for mu = 2,10 and xm = .1,1,2 and κm = 20,200
# (variables are named f1 because it's the first figure of results)

include("larva-base.jl")

μ_vec = [2.,6.,10.]; # Lμ = length(μ_vec);
Lϵ = 20; ϵ_vec = collect(range(1.,0.,length = Lϵ));

# initialize df
df_in = DataFrame(repeat([Float64],3),[:μ0,:ϵ,:μ])
df1 = DataFrame(repeat([Float64],5),[:S,:W,:M,:T,:e])
df2 = DataFrame(repeat([Float64],5),[:S,:W,:M,:T,:e])

for μ0 in μ_vec
    q = copy(p0) 
    q["μ"] = μ0
    M0 = larva_fates(q,other_args).M[1]
    for ϵ in ϵ_vec
        q["ϵ"] = ϵ
        p = copy(q)
        μ = find_zero(
            μ1 -> begin
                pp = copy(p)
                pp["μ"] = μ1
                return(larva_fates(pp,other_args).M[1] - M0)
            end,
            μ0
        )
        p["μ"] = μ
        append!(df_in,DataFrame("μ0" => μ0,"ϵ" => ϵ,"μ" => μ))
        append!(df1,larva_fates(q,other_args))
        append!(df2,larva_fates(p,other_args))
    end
end

rename!(df1,[:S1,:W1,:M1,:T1,:e1])
rename!(df2,[:S2,:W2,:M2,:T2,:e2])
df2[!,:normalized_S2] = df2.S2 ./ (1. .- df2.M2)
X = hcat(df_in,df1,df2)

CSV.write(data * "/S-versus-eps.csv",X)
X = CSV.read(data*"/S-versus-eps.csv")

f2a = plot(
    X,
    x = :ϵ,
    y = :S1,
    color = :μ0,
    Geom.point,
    Geom.line,
    Scale.color_discrete_manual(my_colors...),
    Scale.y_log10,
    Guide.xlabel("Mortality ratio, ε"),
    Guide.ylabel("Prob. of settling, S"),
    Guide.ColorKey("μ"),
    latex_fonts
)
draw(PDF(figz * "/fig-eps-a.pdf",12cm,10cm),f2a)

f2b = plot(
    X,
    x = :ϵ,
    y = :S2,
    color = :μ0,
    Geom.point,
    Geom.line,
    Scale.color_discrete_manual(my_colors...),
    Scale.y_log10,
    Guide.xlabel("Mortality ratio, ε"),
    Guide.ylabel("Prob. of settling, S"),
    latex_fonts
)
draw(PDF(figz * "/fig-eps-b.pdf",12cm,10cm),f2b)


f2c = plot(
    X,
    x = :ϵ,
    y = :μ,
    color = :μ0,
    Geom.point,
    Geom.line,
    Scale.color_discrete_manual(my_colors...),
    Scale.y_log10,
    Guide.xlabel("Mortality ratio, ε"),
    Guide.ylabel("Nearshore mortality, μ"),
    latex_fonts
)
draw(PDF(figz * "/fig-eps-c.pdf",12cm,10cm),f2c)

# hstack(
#     plot(
#         df,
#         x = :ϵ,
#         y = :S1,
#         color = :μ0,
#         Geom.line,
#         Scale.color_discrete_manual(my_colors...),
#         Scale.y_log10
#     ),

#     plot(
#         df,
#         x = :ϵ,
#         y = :S2,
#         color = :μ0,
#         Geom.line,
#         Scale.color_discrete_manual(my_colors...),
#         Scale.y_log10
#     )
# )

# plot(
#     df,
#     x = :ϵ,
#     y = :S,
#     xgroup = :μ0,
#     color = :κm,
#     Geom.subplot_grid(Geom.line,free_y_axis = true),
#     Scale.y_log10
# )

# df1 = filter(row -> (row.χm == 1.),df)

# plot(
#     df1,
#     x = :ϵ,
#     y = :S2,
#     color = :μ0,
#     Scale.color_discrete_manual(my_colors[1:2]...), #(colormap = Colors.colormap("Grays",logscale = true)),
#     Geom.line
# )


# larva_pdf = larva(p0)
# f1_p0 = copy(p0) # base parameter set
# f1_Lϵ = 50 # number of ϵ values to use


# f1_pars = OrderedDict( # parameters to be varied
#     "ϵ" => collect(range(0.,1,length = f1_Lϵ)),
#     "χm" => [.1,1,2],
#     "κm" => [20.,200.],
#     "μ" => [2.,10.]
# )

# # do the parameter sweep, and time it
# f1_time = time()
# f1_df = parsweep_multi(larva_fates,f1_p0,f1_pars,[:S,:D,:W],other_args)
# f1_time = time() - f1_time

# # make figures
# f1a = Gadfly.plot(
#     filter(raw -> (raw.μ < 5),f1_df),
#     y = :S,
#     Scale.y_log10,
#     x = :ϵ,
#     linestyle = :χm,
#     color = :κm,
#     Scale.color_discrete_manual(my_colors[1:2]...), #(colormap = Colors.colormap("Grays",logscale = true)),
#     Geom.line,
#     Coord.cartesian(xmin = 0,xmax = 1,ymin = -3.5,ymax = 0),
#     Guide.xlabel("ε"),
#     Guide.ColorKey(title = "κ")
# )

# f1b = Gadfly.plot(
#     filter(raw -> (raw.μ > 5),f1_df),
#     y = :S,
#     Scale.y_log10,
#     x = :ϵ,
#     linestyle = :χm,
#     color = :κm,
#     Scale.color_discrete_manual(my_colors[1:2]...), #(colormap = Colors.colormap("Grays",logscale = true)),
#     Geom.line,
#     Coord.cartesian(xmin = 0,xmax = 1,ymin = -3.5,ymax = 0),
#     Guide.xlabel("ε"),
#     Guide.ColorKey(title = "κ")
# )
