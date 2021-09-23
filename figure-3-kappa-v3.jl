include("larva-base.jl")

###############################################################
## Figure 3a -- logS v logK for eps = 1 and eps = 0 (three xm's)
a = copy(p0); a["μ"] = 10.;
Lk = 30;
pars_unif = OrderedDict("ϵ" => [1.],"χm" => [1.],"κm" => collect(10. .^range(0.,3,length = Lk)));
pars_hetero = OrderedDict("ϵ" => [0.],"χm" => [.1,1.,2.],"κm" => collect(10. .^range(0.,3,length = Lk)));

A1 = parsweep_multi(larva_fates,a,pars_unif,[:S,:M,:W],other_args);
A2 = parsweep_multi(larva_fates,a,pars_hetero,[:S,:M,:W],other_args);

f3a = plot(
    layer(A1,x = :κm,y = :S,Geom.line,Geom.point,Theme(default_color = my_colors[1])),
    # layer(filter(r -> (r.χm == .1),A2),x = :κm,y = :S,Geom.line,Geom.point,Theme(default_color = my_colors[2])),
    layer(filter(r -> (r.χm == 1.),A2),x = :κm,y = :S,Geom.line,Geom.point,Theme(default_color = my_colors[3])),
    layer(filter(r -> (r.χm == 2.),A2),x = :κm,y = :S,Geom.line,Geom.point,Theme(default_color = my_colors[2])),
    Scale.x_log10,Scale.y_log10,
    latex_fonts,
    Guide.xlabel("Eddy diffusivity, κ̄"),
    Guide.ylabel("Prob. of settling, S")
)
draw(PDF(figz * "/fig-kappa-a-v3.pdf",12cm,10cm),f3a)


###############################################################
## Figures 3b, 3c: eps v K* for xm = 1,2
Le = 25;
Lk = 200;
pars_minmax = OrderedDict("χm" => [1.,2.], "ϵ" => collect(range(0,1,length = Le)),"κm" => collect(10. .^range(0,3,length = Lk)));

df = parsweep_multi(larva_fates,a,pars_minmax,[:S,:M,:W],other_args)

B = DataFrame(repeat([Float64],6),[:ϵ,:χm,:bmax,:bmin,:max,:min])
for χm in pars_minmax["χm"]
    for ϵ in pars_minmax["ϵ"]
        D = filter(row -> ((row.ϵ == ϵ) && (row.χm == χm)),df)
        BMax = ifelse(D.S[1] > D.S[2],1.,NaN)
        BMin = ifelse(D.S[1] < D.S[2],1.,NaN)
        Max = NaN
        Min = NaN
        for j = 2:(Lk-1)
            if D.S[j] < D.S[j-1] && D.S[j] < D.S[j+1]
                Min = pars_minmax["κm"][j]
            end
            if D.S[j] > D.S[j-1] && D.S[j] > D.S[j+1]
                Max = pars_minmax["κm"][j]
            end
        end
        append!(B,[:ϵ => ϵ,:χm => χm,:bmax => BMax,:bmin => BMin,:max => Max,:min => Min])
    end
end
C = stack(B,[:bmin,:bmax,:min,:max],variable_name = :extremum_type,value_name = :extremum)
C[!,:Location] = map(r -> ifelse(string(r)[1] == 'b',:boundary,:interior),C.extremum_type)
C[!,:Type] = map(r -> ifelse(string(r)[(end-2):end] == "min",:min,:max),C.extremum_type)

CSV.write(data*"/critical-kappa-v3.csv",C)

f3b = plot(
    filter(r -> (r.χm == 1.), C)[Int(Lk/2):-1:1,:], # reverse so Max is first
    x = :ϵ,
    y = :extremum,
    shape = :Type,
    color = :Type,
    Geom.point, # Geom.line,
    Scale.y_log10,
    Guide.ylabel("Critical value, κ*"),
    Guide.xlabel("Mortality ratio, ε"),
    Coord.Cartesian(ymax = 3),
    Scale.color_discrete_manual(colorant"black",colorant"white"),
    Theme(discrete_highlight_color=c->colorant"black"),
    latex_fonts2,
)
draw(PDF(figz * "/fig-kappa-b-v3.pdf",13cm,10cm),f3b)


f3c = plot(
    filter(r -> (r.χm == 2.), C)[Int(Lk/2):-1:1,:], # reverse so Max is first,
    x = :ϵ,
    y = :extremum,
    shape = :Type,
    color = :Type,
    Geom.point,
    Scale.y_log10,
    Guide.ylabel("Critical value, κ*"),
    Guide.xlabel("Mortality ratio, ε"),
    Scale.color_discrete_manual(colorant"black",colorant"white"),
    Coord.Cartesian(ymax = 3),
    Scale.color_discrete_manual(colorant"black",colorant"white"),
    Theme(discrete_highlight_color=c->colorant"black"),
    latex_fonts2
)
draw(PDF(figz * "/fig-kappa-c-v3.pdf",13cm,10cm),f3c)
