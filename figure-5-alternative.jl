include("larva-base.jl")

xvec = float.(0:p0["h"]:p0["L"]); NX = length(xvec);
tvec = float.(0:p0["k"]:1); NT = length(tvec);

p_wu = copy(p0); p_wu["μ"] = 2.;
p_wh = copy(p0); p_wh["μ"] = 2.; p_wh["ϵ"] = 0.;
p_su = copy(p0); p_wu["μ"] = 10.;
p_sh = copy(p0); p_wh["μ"] = 10.; p_wh["ϵ"] = 0.;

Q_wu,_ = larva(p_wu,other_args...);
Q_wh,_ = larva(p_wh,other_args...);
Q_su,_ = larva(p_su,other_args...);
Q_sh,_ = larva(p_sh,other_args...);

Q_wu = DataFrame(:x => repeat(xvec,outer = NT),:t => repeat(tvec,inner = NX),:p => vec(Q_wu));
Q_wh = DataFrame(:x => repeat(xvec,outer = NT),:t => repeat(tvec,inner = NX),:p => vec(Q_wh));
Q_su = DataFrame(:x => repeat(xvec,outer = NT),:t => repeat(tvec,inner = NX),:p => vec(Q_su));
Q_sh = DataFrame(:x => repeat(xvec,outer = NT),:t => repeat(tvec,inner = NX),:p => vec(Q_sh));

set_default_plot_size(20cm,20cm)
hstack(
    vstack(
        plot(Q_wu,x = :t,y = :x,color = :p,Geom.rectbin,Coord.cartesian(xmin = 0,xmax = 1,ymin = 0,ymax = 10)),
        plot(Q_su,x = :t,y = :x,color = :p,Geom.rectbin,Coord.cartesian(xmin = 0,xmax = 1,ymin = 0,ymax = 10))
    ),
    vstack(
        plot(Q_wh,x = :t,y = :x,color = :p,Geom.rectbin,Coord.cartesian(xmin = 0,xmax = 1,ymin = 0,ymax = 10)),
        plot(Q_sh,x = :t,y = :x,color = :p,Geom.rectbin,Coord.cartesian(xmin = 0,xmax = 1,ymin = 0,ymax = 10))
    )
)
