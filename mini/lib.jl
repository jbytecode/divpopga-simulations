import DivPopGa.ClusteredGa as CLGA
using DataFrames
using CSV
using Statistics
using HypothesisTests
using RCall
include("fun.jl")

function functionLoader(funName)
    return Symbol("fun_",funName)
end

function mannwhitney(x, y; alternative = "greater")
    @rput x
    @rput y
    if !(alternative in ["less", "greater", "two.sided"])
        @error "alternative should be one of the 'less', 'greater', or 'two.sided'"
        return -1 
    end
    result = R"wilcox.test(x, y, paired = FALSE, alternative = $(alternative))"
    pvalue = result[3]
    return convert(Float64, pvalue)
end