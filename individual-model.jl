# THIS FILE CONTAINS CODE FOR THE SDE MODEL OF AN INDIVIDUAL LARVA.

# functions for 1 time-step using Euler-Maruyama and stochastic Runge-Kutta
# p = parameter dictionary; f = drift; g = √diffusivity; x = current location;
# t = current time; dt = time-step length; dB = brownian increment over dt
em1d_step(p,f,g,x,t,dt,dB) = x + f(t,x,p)*dt + g(t,x,p)*dB

function srk1d_step(p,f,g,x,t,dt,dB)
    F = f(t,x,p)
    G = g(t,x,p)
    dt2 = sqrt(dt)
    Z = x + F*dt + G*dt2
    return(x + F*dt + G*dB + (g(t,Z,p) - G)*(dB^2 - dt)/(2*dt2))
end

# simulates 1d killed diffusion process starting at x0 on time mesh tvec.
# p = parameters; f,g = drift and √diffusivity; x0 = initial value;
# tvec = time mesh; BCs = dictionary specifying boundary conditions, of form
#   BCs[:left or :right] = [value,type] where type = :a (absorbing) or :r (reflecting)
# k = killing rates; kill_names = vector of names for each type of killing;
# step_method = function for simulating 1 time step (defaults to Euler-Maruyama)
function KilledDiffusion1d(p,f,g,x0,tvec,BCs,k,kill_names;step_method = em1d_step)
    # initialize
    N = length(tvec)
    X = x0*ones(N)
    flag = "ran to tmax"
    # compute random numbers
    dB = randn(N-1) # brownian increments
    U = rand(N-1) # uniform rv's to decide which killing (if any) occurs

    # do the simulation
    for n = 1:(N-1)
        x = X[n]
        t = tvec[n]
        dt = tvec[n+1] - t
        kill_rates = k(t,x,p)

        # check if survived or killed
        if U[n] > 1 - exp(-dt*sum(kill_rates))
            # if survived, compute next step...
            Z = step_method(p,f,g,x,t,dt,dB[n])

            # ... and see if we've hit a boundary
            if Z < BCs[:left][1]
                if BCs[:left][2] == :a
                    X[(n+1):N] = BCs[:left][1]
                    flag = "absorbed on left"
                    break
                else
                    X[n+1] = BCs[:left][1] + abs(Z - BCs[:left][1])
                end
            elseif Z > BCs[:right][1]
                if BCs[:right][2] == :a
                    X[(n+1):N] = BCs[:right][1]
                    flag = "absorbed on right"
                    break
                else
                    X[n+1] = BCs[:right][1] - abs(Z - BCs[:right][1])
                end
            else
                X[n+1] = Z
            end
        else # if killed, decide what killed it and end process
            # kill_probs = kill_rates/sum(kill_rates)
            flag = sample(kill_names,Weights(kill_rates))
            X[(n+1):N] .= NaN
            break
        end
    end
    return(X,flag)
end

# runs M copies of KilledDiffusion1d and reports outcomes. If kill_type specified,
# runs until M trials end in specified type
function KilledDiffusion1d_many(M,p,f,g,X0,tvec,BCs,k,kill_names;step_method = em1d_step,kill_type = nothing)
    t0 = time()
    N = length(tvec)
    XX = zeros(N,M)
    flags = Vector{String}(undef,M)
    if length(X0) == 1
        X0 = ones*X0
    end

    if isnothing(kill_type) == true
        # if no kill_type specified, run in parallel (exactly M sims)
        Threads.@threads for m = 1:M
            XX[:,m],flags[m] = KilledDiffusion1d(p,f,g,X0[m],tvec,BCs,k,kill_names;step_method = step_method)
        end
    else
        # if kill_type given, run until M of those happen (not parallelized)
        m = 0
        while m < M
            X,flag = KilledDiffusion1d(p,f,g,X0[m+1],tvec,BCs,k,kill_names;step_method = step_method)
            if flag == kill_type
                m += 1
                XX[:,m],flags[m] = X,flag
            end
        end
    end

    return(XX,flags,-t0 + time())
end

# compute path stats (distance traveled, max value attained, time spent with X > 1)
# based on M killed diffusions
function KilledDiffusion1d_DistTraveled(M,p,f,g,x0,tvec,BCs,k,kill_names;step_method = em1d_step,kill_type = nothing,time_beyond_D = [1])
    dt = tvec[2] - tvec[1]
    XX,_,Time = KilledDiffusion1d_many(M,p,f,g,x0,tvec,BCs,k,kill_names;step_method = step_method,kill_type = kill_type)
    df = DataFrame(
        :TotalDistance => [sum(filter(!isnan,XX[:,m]))*dt for m = 1:M],
        :MaxDistance => [maximum(filter(!isnan,XX[:,m])) for m = 1:M],
        # :TimeBeyondHab => [sum(filter(!isnan,XX[:,m]) .> 1)*dt for m = 1:M]
    )
    for D in time_beyond_D
        df[!,Symbol("TimeBeyond"*string(Int(round(D))))] = [sum(filter(!isnan,XX[:,m]) .> D)*dt for m = 1:M]
    end
    return(XX,df,Time)
end
