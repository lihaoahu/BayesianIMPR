module BayesianIPR

using Turing
using Distributions
using LinearAlgebra
using StatsPlots

export 
    # modules
    AHP,
    IntervalAHP,

    # types
    IMPR,
    IAPR,
    MPR,
    Interval,

    # functions
    getdata,
    getall,
    reorder,
    ahp,
    ahp_analy,
    iahp,
    iahp_analy,
    iahp_analy_uninfo,
    acceptability_index_mipr




include("IntervalAnalysis.jl")
include("PreferenceRelations.jl")
include("AHP.jl")
include("IntervalAHP.jl")



end # module BayesianIPR
