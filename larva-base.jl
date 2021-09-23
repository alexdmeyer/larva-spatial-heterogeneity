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
using Roots
using SpecialFunctions
using CSV

home = "/Users/alexandermeyer/Dropbox/Research/Project01_larvalSuccess/Manuscript/Revisions-2021-01-11/larval-delivery-revisions"
figz = home * "/figures-revised"
data = home * "/data"

include("pop-model.jl")
include("individual-model.jl")
include("parameter-sweep.jl")
include("bisection.jl")

# Movement functions and parameters
κ(τ,χ,p) = p["κm"] * ((χ/p["χb"]).^p["α"] .* (χ .< p["χb"]) .+ (χ .>= p["χb"]))
κχ(τ,χ,p) = p["κm"] * p["α"]*(χ/p["χb"]).^(p["α"]-1) .* (χ .< p["χb"])
ν(τ,χ,p) = zeros(size(χ));
movepars = Dict("κm" => 100., "χb" => 0., "α" => 2.)

# settling function and parameters
λs(τ,χ,p) = p["σ"] * (τ .>= p["τc"]) .* (χ .<= 1);
settlepars = Dict("σ" => 2., "τc" => .5)

# mortality functions and parameters
λd(τ,χ,p) = p["μ"] * (p["ϵ"] .+ (1 - p["ϵ"])*(χ .<= p["χm"]))
mortpars = Dict("μ" => 6. ,"ϵ" => 1.,"χm" => 1.)

# additional parameters
δ = .1 # basis for time and space meshes
otherpars = Dict("L" => 2*movepars["κm"]^.5*erfinv(.99), "h" => .5δ^2, "k" => δ,"χ0" => nothing)

# parameter dict, ranges, additional arguments
other_args = (κ,ν,λd,λs,"d","d") #(κ,ν,λd,λs)
p0 = merge(movepars,settlepars,mortpars,otherpars)

# for stochastic DE model
BCs = Dict(:left => [0,:r],:right => [Inf,:r])
f(t,x,p) = ν(t,x,p)[1] + κχ(t,x,p)[1]
g(t,x,p) = sqrt(2*κ(t,x,p)[1])
k(t,x,p) = [λd(t,x,p), λs(t,x,p)]
kill_names = ["death","settling"]

# stuff for plotting
gadfly_colors = Scale.color_discrete_hue().f(10)
my_color_subset = [1 6 4 5]
my_colors = gadfly_colors[my_color_subset]
ρ = 0.025 # small radius for searching dataframes for parameter values
# (don't want to use μ == ? if there's even a little floating point error)
latex_fonts = Theme(major_label_font="CMU Serif", major_label_font_size=16pt,
                    minor_label_font="CMU Serif", minor_label_font_size=14pt,
                    key_title_font="CMU Serif", key_title_font_size=12pt,
                    key_label_font="CMU Serif", key_label_font_size=10pt);
latex_fonts2 = Theme(major_label_font="CMU Serif", major_label_font_size=16pt,
                    minor_label_font="CMU Serif", minor_label_font_size=14pt,
                    key_title_font="CMU Serif", key_title_font_size=12pt,
                    key_label_font="CMU Serif", key_label_font_size=10pt,
                    discrete_highlight_color=c->colorant"black");