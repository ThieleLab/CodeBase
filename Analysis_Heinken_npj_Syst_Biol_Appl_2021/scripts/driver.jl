#-------------------------------------------------------------------------------------------
#=
    Purpose:    driver to run a set of models
    Author:     Laurent Heirendt - LCSB - Luxembourg
    Date:       December 2016
=#

#-------------------------------------------------------------------------------------------

# Analysis performed
 compute all exchanges and microbe contributions to overall metabolite production 
 Loops through all internal exchanges (time-consuming)

# set the directory for storing results
    resultsDir = ENV["HOME"]*"/results_MicrobeContributions"  # "/

# create the directory if it does not exist
if !isdir(resultsDir) mkdir(resultsDir) end

# define model and solution parameters
    optPercentage = 0.0

# set the name of the stoichiometric matrix
    maxtrixAS = "A"

# set the solver name
solverName        = :CPLEX

# set solver parameters
    solParams = [
        # decides whether or not results are displayed on screen in an application of the C API.
        (:CPX_PARAM_SCRIND,         0);

        # sets the parallel optimization mode. Possible modes are automatic, deterministic, and opportunistic.
        (:CPX_PARAM_PARALLELMODE,   1);

        # sets the default maximal number of parallel threads that will be invoked by any CPLEX parallel optimizer.
        (:CPX_PARAM_THREADS,        1);

        # partitions the number of threads for CPLEX to use for auxiliary tasks while it solves the root node of a problem.
        (:CPX_PARAM_AUXROOTTHREADS, 2);

        # decides how to scale the problem matrix.
        (:CPX_PARAM_SCAIND,         -1);

        # controls which algorithm CPLEX uses to solve continuous models (LPs).
        (:CPX_PARAM_LPMETHOD,       0)
    ] #end of solParams

# set the model name
    modelName         = "model"

#set the benchmark boolean
benchmark         = true

# define the splitting strategy
strategy          = 0

# set the number of workers on host node (local node)
    nWorkersLoc = 24

# include as well SSH nodes and launch remote processes
connectSSHWorkers = false # true | false

#-------------------------------------------------------------------------------------------
# load the models from the directory on the server

# set the local path for public cobratoolbox and model
    MODEL_PATH = ENV["HOME"]*"/$(gethostname())/"

# read the directory with all the models
modelDir = readdir(MODEL_PATH)

# retrieve the number of models
nModels = length(modelDir)
info("The analysis will be performed on $nModels model(s) ...\n")

#-------------------------------------------------------------------------------------------

# load the connection details
include("$(Pkg.dir("COBRA"))/src/connect.jl")

# create a parallel pool and determine its size
if isdefined(:nWorkersLoc) && isdefined(:connectSSHWorkers)
    workersPool, nWorkers = createPool(nWorkersLoc, connectSSHWorkers, "sshCfgPrivate.jl")
end

@everywhere using MAT
@everywhere using COBRA

# set the COBRA solver
solver = changeCobraSolver(solverName, solParams);

# intialize a vector to store all benchmark times
bmTimes = zeros(length(modelDir))

# loop through all the models
for i = 1:nModels

    global minFlux, maxFlux, optSol, fbaSol, fvamin, fvamax, statussolmin, statussolmax, rxnsList, rxnsOptMode, sumResults

    info("Loading $MODEL_PATH/$(modelDir[i]) ...\n")

    # load the model i
	model = loadModel("$MODEL_PATH$(modelDir[i])", maxtrixAS, modelName)

    # define an empty reaction list vector
        rxnsList = []

        # find which reactions are internal exchange reactions for the target metabolite
        counter = 0
        pattern = "Diet_EX_"
        for j = 1:length(model.rxns)
            if contains(model.rxns[j],pattern)
             push!(rxnsList, j)
                counter += 1
            end
        end
        pattern = "_IEX_"
        for j = 1:length(model.rxns)
            if contains(model.rxns[j],pattern)
             push!(rxnsList, j)
                counter += 1
            end
        end
        pattern = "EX_"
        for j = 1:length(model.rxns)
            if contains(model.rxns[j],pattern)
             push!(rxnsList, j)
                counter += 1
            end
        end


   # initialize a vector to store all minFlux and maxFlux values

        sumResults = zeros(counter, 5)
        info("sumResults initialized")

    # save chunks of the fvamin and fvamax vectors
    saveChunks = true

    # set the optimization mode for each reaction to 2 (both maximization and minimization)

        rxnsOptMode = 2 + zeros(Int64, length(rxnsList))
        info("There are $counter reaction(s) that will be maximized and minimized.\n")

    # launch the distributedFBA for the given list of reactions
    startTime   = time()

        minFlux, maxFlux, optSol, fbaSol, fvamin, fvamax, statussolmin, statussolmax = distributedFBA(model, solver, nWorkers=nWorkers, optPercentage=optPercentage, objective=objective, rxnsList=rxnsList, strategy=strategy, rxnsOptMode=rxnsOptMode, preFBA=false, saveChunks=saveChunks, resultsDir=resultsDir, logFiles=false, onlyFluxes=true)

    solTime = time() - startTime

    # output and save the flux values

        sumResults[:, 1] = rxnsList
        sumResults[:, 2] = minFlux[rxnsList]
        sumResults[:, 3] = maxFlux[rxnsList]
        sumResults[:, 4] = statussolmin[rxnsList]
        sumResults[:, 5] = statussolmax[rxnsList]

    # output solution summary
    printSolSummary(modelDir[i][1:end-4], optSol, maxFlux, minFlux, solTime, nWorkers, solverName, strategy, saveChunks)

    # save all the calculated variables as a .mat file

        saveDistributedFBA("$(resultsDir)/summary_$(modelDir[i][1:end-4]).mat", ["sumResults", "rxnsList", "minFlux", "maxFlux", "statussolmin", "statussolmax"])
    
    end

    # report and save the benchmark solution time
    if benchmark
        bmTimes[i] = solTime
        print_with_color(:yellow, "Solution time for benchmark $(bmTimes[i])\n\n")
    end
end

# print out a summary of the benchmark times for all models
print_with_color(:red, "--- Benchmark times with $nWorkers workers ---\n ")
for i = 1:nModels
    print_with_color(:cyan, "$(modelDir[i][1:end-4])")
    print("       $(bmTimes[i]) [s]")
    print("\n")
end
