# General function for computing outputs of interest (in this case S,D,W)
# from Model (which takes arguments p, other_args; in this case, larva_fates)
# over an (ordered) dictionary of values for a few focal parameters.

# Model = function for outputs; p0 = base parameter values; pDict = values of
# parameters to vary; outputs = subset of outputs we actually care about;
# other_args = arguments besides parameters p for Model

# output is a dataframe whose first Npars columns are parameter values
# used (ignoring unvaried ones) and last few columns are outputs of interest
function parsweep_multi(Model,p0,pDict,outputs,other_args)
    # 1. create dataframe of parameter combinations
    parNames = collect(keys(pDict))
    Nvals_each = [length(pDict[k]) for k in parNames]
    Npars = length(Nvals_each)
    Nvals = prod(Nvals_each)
    Inner = cumprod([1;Nvals_each[1:(Npars-1)]])
    Outer = Int.(Nvals./cumprod(Nvals_each))

    df_pars = DataFrame()
    for n in 1:Npars
        parStr = parNames[n]
        parSym = Symbol(parStr)
        df_pars[!,parSym] = repeat(pDict[parStr],inner = Inner[n],outer = Outer[n])
    end

    # 2. for each row of the dataframe, update the par dictionary
    # and compute outputs, writing those into another dataframe
    df_outs = DataFrame(OrderedDict(oname => [] for oname in outputs))
    for n = 1:nrow(df_pars)
        this_p = copy(p0)
        for par in names(df_pars)
            this_p[String(par)] = df_pars[n,par]
        end
        # print(this_p)
        model_out = Model(this_p,other_args)
        append!(df_outs,model_out[:,outputs])
    end

    # 3. merge the two df's into a single output
    return(hcat(df_pars,df_outs))
end
