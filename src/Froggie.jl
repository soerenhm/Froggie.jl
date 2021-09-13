module Froggie

using AxisArrays, Unitful
const axes = Base.axes  # due to possible conflict with AxisArrays.axes

include("utils.jl")
include("core.jl")
include("binning.jl")

export frogtrace, freqdomain, binalong

end # module
