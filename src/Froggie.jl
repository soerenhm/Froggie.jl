module Froggie

using AxisArrays, Unitful
const axes = Base.axes  # due to possible conflict with AxisArrays.axes


include("binning.jl")

export binalong

end # module
