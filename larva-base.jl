# THIS FILE CREATES PARAMETER DICTIONARIES AND FUNCTIONS FOR THE LARVA MODEL.
# IT ALSO IMPORTS ALL THE PACKAGES WE WILL NEED AND RUNS THE FILES WHICH
# CREATE THE ACTUAL MODEL FUNCTIONS

using Pkg
using SparseArrays
using DataFrames
using DataStructures
using Random
using StatsBase
using Distributions
using Gadfly, Cairo, Fontconfig
include("pop-model.jl")
include("individual-model.jl")
include("parameter-sweep.jl")

# Movement functions and parameters
κ(t,x,p) = p["κ̄"] * ((x/p["xb"]).^p["α"] .* (x .< p["xb"]) .+ (x .>= p["xb"]))
κx(t,x,p) = p["κ̄"] * p["α"]*(x/p["xb"]).^(p["α"]-1) .* (x .< p["xb"])
ν(t,x,p) = zeros(size(x));
movepars = Dict("κ̄" => 100., "xb" => 0., "α" => 2.)

# settling function and parameters
λs(t,x,p) = p["σ"] * (t .>= p["tpc"]) .* (x .<= 1);
settlepars = Dict("σ" => 2., "tpc" => .5)

# mortality functions and parameters
λm(t,x,p) = p["μ"] * (p["ϵ"] .+ (1 - p["ϵ"])*(x .<= p["xm"]))
mortpars = Dict("μ" => 6. ,"ϵ" => 1.,"xm" => 1.)

# additional parameters
δ = .1 # basis for time and space meshes
otherpars = Dict("L" => 30., "h" => .5δ^2, "k" => δ,"x0" => nothing)

# parameter dict, ranges, additional arguments
other_args = (κ,ν,λm,λs,"d","d") #(κ,ν,λd,λs)
p0 = merge(movepars,settlepars,mortpars,otherpars)

# for stochastic DE model
BCs = Dict(:left => [0,:r],:right => [Inf,:r])
f(t,x,p) = ν(t,x,p)[1] + κx(t,x,p)[1]
g(t,x,p) = sqrt(2*κ(t,x,p)[1])
k(t,x,p) = [λm(t,x,p), λs(t,x,p)]
kill_names = ["mortality","settling"]

# stuff for plotting
gadfly_colors = Scale.color_discrete_hue().f(10)
my_color_subset = [1 6 4 5]
my_colors = gadfly_colors[my_color_subset]
ρ = 0.025 # small radius for searching dataframes for parameter values
# (don't want to use μ == ? if there's even a little floating point error)
