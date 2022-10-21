
import Metaheuristics
using Optim 
import SQLite
import SQLite.DBInterface
import Base.Threads

include("../lib/test-suite.jl")


# Params 
popsize = 100
p = 10
lowers = repeat([-100.0], outer = p)
uppers = repeat([100.0], outer = p)

# THE FUNCTIONS
costfn = cec13
costfnname = "cec13"
# END


# DB 
db = SQLite.DB("sim.db")
SQLite.execute(db, 
    """
    CREATE TABLE IF NOT EXISTS Simulation 
    (fname text, p integer, de float, pso float, abc float, eca float, woa float, mccga float)
    """
)

for sims = 1:100
    # DE 
    res_de = Threads.@spawn Metaheuristics.optimize(costfn, [lowers uppers], Metaheuristics.DE(N = popsize))
    # PSO 
    res_pso = Threads.@spawn Metaheuristics.optimize(costfn, [lowers uppers], Metaheuristics.PSO(N = popsize))
    # ABC 
    res_abc = Threads.@spawn Metaheuristics.optimize(costfn, [lowers uppers], Metaheuristics.ABC(N = popsize))
    # ECA 
    res_eca = Threads.@spawn Metaheuristics.optimize(costfn, [lowers uppers], Metaheuristics.ECA(N = popsize))
    # WOA
    res_woa = Threads.@spawn Metaheuristics.optimize(costfn, [lowers uppers], Metaheuristics.WOA(N = popsize))
    # MCCGA 
    res_mccga = Threads.@spawn Metaheuristics.optimize(costfn, [lowers uppers], Metaheuristics.MCCGA(N = popsize))

    de    = fetch(res_de).best_sol.f
    pso   = fetch(res_pso).best_sol.f
    abc   = fetch(res_abc).best_sol.sol.f
    eca   = fetch(res_eca).best_sol.f
    woa   = fetch(res_woa).best_sol.f
    mccga = fetch(res_mccga).best_sol.f

    sql = """
    INSERT INTO Simulation 
    (fname, p, de, pso, abc, eca, woa, mccga) 
    VALUES
    ('$costfnname', $p, $de, $pso, $abc, $eca, $woa, $mccga)
    """
    SQLite.DBInterface.execute(db, sql)
    printstyled(".", color = :green)
end 

println()

SQLite.DBInterface.close!(db)