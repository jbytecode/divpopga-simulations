using Pkg
Pkg.develop(path = "../divpopga/.")
import DivPopGa.ClusteredGa as CLGA
import Metaheuristics
using Optim 
using Revise

# BENT CIGAR FUNCTION
function cec1(x::Array{Float64, 1})::Float64
    
    term1 = x[1]^2
    term2 = 0
    for i in 2:length(x)
        term2 = term2 + x[i]^2
    end
    term2 = term2 * 1000000

    return term1 + term2

end

# Params 
popsize = 200
generations = 5000
mutationParameter1 = 3.0
mutationParameter2 = 0.01
num_elitists = 2
lowerbound = -100.0
upperbound = 100.0
numV = 10

# THE FUNCTIONS
costfn = cec1
crossfn = CLGA.makeblxalphacrossover()
mutatefn  = CLGA.makenormalmutation(mutationParameter1, mutationParameter2)
# END

# CLASSIC GA
res_ga = CLGA.ga(
    popsize,
    generations,
    repeat([lowerbound], outer=numV),      
    repeat([upperbound], outer=numV),      
    costfn,     
    crossfn,    
    mutatefn,   
    CLGA.GA_TYPE_CLASSIC, 
    elitism = num_elitists
)

# KMEANS
res_kmeans = CLGA.ga(
    popsize,    
    generations,
    repeat([lowerbound], outer=numV),      
    repeat([upperbound], outer=numV),      
    costfn,     
    crossfn,    
    mutatefn,   
    CLGA.GA_TYPE_CLUSTER, 
    elitism = num_elitists
)

# KMEANSSIM
res_kmeanssim = CLGA.ga(
    popsize,    
    generations,
    repeat([lowerbound], outer=numV),      
    repeat([upperbound], outer=numV),      
    costfn,     
    crossfn,    
    mutatefn,   
    CLGA.GA_TYPE_CLUSTER_SIM, 
    elitism = num_elitists
)

# DE 
res_de = Metaheuristics.optimize(costfn, [lowers uppers], DE())

# PSO 
res_pso = Metaheuristics.optimize(costfn, [lowers uppers], PSO())

# MCCGA 
res_mccga = Metaheuristics.optimize(costfn, [lowers uppers], MCCGA())
