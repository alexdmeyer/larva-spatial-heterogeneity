include("larva-base.jl")

LX = 20;
Case1 = OrderedDict( # moderate mortality, uniform in space
    "μ" => [6.],
    "ϵ" => [1.],
    "χm" => [1.],
    "χb" => vcat(0,collect(range(.75,10,length = LX-1)))
) 
df1 = parsweep_multi(larva_fates,p0,Case1,[:S,:M,:W,:e],other_args)

Case2 = OrderedDict(# strong mortality, medium HMZ
    "μ" => [6.],
    "ϵ" => [0.],
    "χm" => [1.],
    "χb" => vcat(0,collect(range(.75,10,length = LX-1)))
)
df2 = parsweep_multi(larva_fates,p0,Case2,[:S,:M,:W,:e],other_args)

Case3 = OrderedDict(# strong mortality, medium HMZ
    "μ" => [6.],   
    "ϵ" => [0.],
    "χm" => [3.],
    "χb" => vcat(0,collect(range(.75,10,length = LX-1)))
)
df3 = parsweep_multi(larva_fates,p0,Case3,[:S,:M,:W,:e],other_args)

Y = vcat(df1,df2,df3);
CSV.write(data * "/CBL-Data.csv",Y)



# plot(df,x = :χb,y = :S,color = :χm,linestyle = :ϵ,Geom.line,Geom.point,latex_fonts,Guide.xlabel("CBL width, xb"),Guide.ylabel("Prob. of settling, S"))
# plot(df,x = :χb,y = :M,color = :χm,linestyle = :ϵ,Geom.line,Geom.point,latex_fonts,Guide.xlabel("CBL width, xb"),Guide.ylabel("Prob. of mortality, M"))
# plot(df,x = :χb,y = :W,color = :χm,linestyle = :ϵ,Geom.line,Geom.point,latex_fonts,Guide.xlabel("CBL width, xb"),Guide.ylabel("Prob. of wastage, W"))

f4a = plot(
    layer(filter(r -> (r.ϵ == 1.),Y),x = :χb,y = :S,Geom.line,Geom.point,Theme(default_color = my_colors[1])),
    layer(filter(r -> (r.ϵ == 0. && r.χm == 1.),Y),x = :χb,y = :S,Geom.line,Geom.point,Theme(default_color = my_colors[3])),
    layer(filter(r -> (r.ϵ == 0. && r.χm == 3.),Y),x = :χb,y = :S,Geom.line,Geom.point,Theme(default_color = my_colors[2])),
    latex_fonts,
    Guide.xlabel("CBL width, xb"),
    Guide.ylabel("Prob. of settling, S")
)
draw(PDF(figz * "/fig-cbl-a.pdf",12cm,10cm),f4a)


f4b = plot(
    layer(filter(r -> (r.ϵ == 1.),Y),x = :χb,y = :M,Geom.line,Geom.point,Theme(default_color = my_colors[1])),
    layer(filter(r -> (r.ϵ == 0. && r.χm == 1.),Y),x = :χb,y = :M,Geom.line,Geom.point,Theme(default_color = my_colors[3])),
    layer(filter(r -> (r.ϵ == 0. && r.χm == 3.),Y),x = :χb,y = :M,Geom.line,Geom.point,Theme(default_color = my_colors[2])),
    latex_fonts,
    Guide.xlabel("CBL width, xb"),
    Guide.ylabel("Prob. of mortality, M")
)
draw(PDF(figz * "/fig-cbl-b.pdf",12cm,10cm),f4b)


f4c = plot(
    layer(filter(r -> (r.ϵ == 1.),Y),x = :χb,y = :W,Geom.line,Geom.point,Theme(default_color = my_colors[1])),
    layer(filter(r -> (r.ϵ == 0. && r.χm == 1.),Y),x = :χb,y = :W,Geom.line,Geom.point,Theme(default_color = my_colors[3])),
    layer(filter(r -> (r.ϵ == 0. && r.χm == 3.),Y),x = :χb,y = :W,Geom.line,Geom.point,Theme(default_color = my_colors[2])),
    latex_fonts,
    Guide.xlabel("CBL width, xb"),
    Guide.ylabel("Prob. of wastage, W")
)
draw(PDF(figz * "/fig-cbl-c.pdf",12cm,10cm),f4c)

