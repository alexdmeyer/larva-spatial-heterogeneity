# THIS FILE CONTAINS THE PDE MODEL FOR THE PROBABILITY DENSITY FUNCTION
# OF LARVA POSITION OVER TIME (AND A FUNCTION FOR JUST COMPUTING PROBABILITIES
# OF FATES, SINCE WE OFTEN DON'T NEED THE FULL DENSITY)

# this is the actual model
# p = parameter dictionary; κ = eddy diffusivity function; ν = advection
# function (usually 0); λm = mortality rate function; λs = settling rate function;
# lbc and rbc = left and right boundary types ("n" = neumann/reflecting,
# "d" = dirichlet/absorbing). Code automatically corrects x = 0 (shoreline)
# boundary conditions if it makes no sense (larvae are lost or division by 0)
function larva(p,κ,ν,λm,λs,lbc,rbc)
    # start timer
    t0 = time()

    # create meshes & extract useful quantities
    xmesh = collect(0:p["h"]:p["L"])
    NX = length(xmesh) - 2
    x_mp = xmesh[2:end] .- p["h"]/2 # midpoints of spatial mesh, χ_mp[n] = h(n-.5)
    tmesh = collect(0:p["k"]:1) # non-dimensionalization hard-coded here
    NT = length(tmesh)

    # initialize -- larvae spawned uniformly over habitat, or concentrated at χ0
    if isnothing(p["x0"])
        # uniform IC
        Q_inner = float.(xmesh[2:(end-1)] .<= 1)
    else
        # point-mass IC
        Q_inner = (abs.(xmesh[2:(end-1)] .- p["x0"]) .< p["h"]/2) / p["h"] # solution on inner mesh at current t
    end

    # pre-allocate
    Q = zeros(NX+2,NT) # full solution, including boundarys
    dead_jd = zeros(NX+2,NT) # joint density of death time/location
    settled_jd = zeros(NX+2,NT) # joint density of settling time/location
    # julia really has no "eye" and "speye"??
    I = sparse([1. * (i == j) for i = 1:NX,j = 1:NX])
    # A is the matrix for our discretized linear PDE. Preallocate as sparse I
    A = copy(I)

    # pre-compute function values on meshes
    κ_mesh = κ(p["k"],xmesh,p)
    κ_mp = κ(p["k"],x_mp,p)
    ν_mesh = ν(p["k"],xmesh,p)
    ν_mp = ν(p["k"],x_mp,p)
    λm_mesh = λm(p["k"],xmesh,p)
    λs_mesh = λs(p["k"],xmesh,p)

    # correct left boundary condition, if needed. Uncomment "println"s if you
    # want Julia to tell you when it does this (but it'll slow down code)
    if lbc == "n" && p["h"]*ν_mesh[1] + κ_mesh[1] < 10^-10
        lbc = "d"
        # println("LBC: no-flux boundary is inappropriate because flux vanishes at boundary anyway. Switching to absorbing boundary.")
    elseif lbc == "d" && p["h"]*ν_mesh[1] + κ_mesh[1] > 10^-10
        lbc = "n"
        # println("LBC: absorbing boundary is inappropriate because larvae are lost at shore. Switching to reflecting boundary.")
    end

    # initialize, then do the solve
    Q[2:(NX+1),1] = Q_inner
    for a in 2:NT
        # compute the next step
        A_update!(A,NX,p["h"],κ_mesh,κ_mp,ν_mesh,ν_mp,λm_mesh,λs_mesh,lbc,rbc)
        Q_inner = (I - p["k"]*A) \ Q_inner

        # store solution and compute next steps for JDs
        Q[2:(NX+1), a] = Q_inner
        # add reflecting boundary condition, if applicable
        if lbc == "n"
            Q[1,a] = Q[2,a] * κ_mesh[1]/(κ_mesh[1] + p["h"]*ν_mesh[1])
        end
        if rbc == "n"
            Q[NX+2,a] = Q[NX+1,a] * κ_mesh[NX+2]/(κ_mesh[NX+2] - p["h"]*ν_mesh[NX+2])
        end

        # update joind densities
        dead_jd[:, a] = Q[:,a] .* λm_mesh
        settled_jd[:, a] = Q[:,a] .* λs_mesh

        # update pre-computed function vectors
        if a < NT
            κ_mesh = κ(a*p["k"],xmesh,p)
            κ_mp = κ(a*p["k"],x_mp,p)
            ν_mesh = ν(a*p["k"],xmesh,p)
            ν_mp = ν(a*p["k"],x_mp,p)
            λm_mesh = λm(a*p["k"],xmesh,p)
            λs_mesh = λs(a*p["k"],xmesh,p)
        end
    end
    # dictionary of joint densities
    jd = Dict("planktonic" => Q, "dead" => dead_jd, "settled" => settled_jd)
    # prop in each class over time (planktonic, dead, settled)
    # Final values are S, M, W
    fxnl = Dict(
        "planktonic" => reshape(sum(Q,dims = 1),NT)*p["h"],
        "dead" => cumsum(reshape(sum(dead_jd,dims = 1),NT))*p["h"]*p["k"],
        "settled" => cumsum(reshape(sum(settled_jd,dims = 1),NT))*p["h"]*p["k"]
    )
    t1 = time() - t0
    return(Q,tmesh,xmesh,jd,fxnl,t1)
end

# function for updating the matrix A each step
function A_update!(A,NX,h,κ_mesh,κ_mp,ν_mesh,ν_mp,λm_mesh,λs_mesh,lbc,rbc)
    # "mp" maps n -> (n - .5)h; "mesh" maps n -> (n-1)h
    for n in 1:(NX-1)
        A[n,n] = -(κ_mp[n] + κ_mp[n+1])/h^2 - λm_mesh[n+1] - λs_mesh[n+1]
        A[n+1,n] = ((ν_mp[n+1] - .5*ν_mesh[n+1]) + κ_mp[n+1]/h)/h
        A[n,n+1] = (-(ν_mp[n+1] - .5*ν_mesh[n+2]) + κ_mp[n+1]/h)/h
    end
    A[NX,NX] = -(κ_mp[NX] + κ_mp[NX+1])/h^2 - λm_mesh[NX+1] - λs_mesh[NX+1]
    if lbc == "n"
        A[1,1] += (κ_mp[1]/h^2 - .5*(2*ν_mp[1] - ν_mesh[2])/h) * (κ_mesh[1]/(h*ν_mesh[1] + κ_mesh[1]))
    end
    if rbc == "n"
        A[NX,NX] += (κ_mp[NX+1]/h^2 + .5*(2*ν_mp[NX+1] - ν_mesh[NX+1])/h) * (κ_mesh[NX+2]/(κ_mesh[NX+2] - h*ν_mesh[NX+2]))
    end
end

# just computes S, W, D. other_args = (κ,ν,λd,λs,lbc,rbc)
function larva_fates(p,other_args)
    _,_,_,_,fx,_ = larva(p,other_args...)
    Y = DataFrame(
        S = fx["settled"][end], # settled
        W = fx["planktonic"][end], # wasted
        M = fx["dead"][end] # dead
    )
    Y.T = Y.W + Y.S + Y.M
    Y.e = Y.T .- 1.
    return(Y)
end
