module Froggie

using Reexport

@reexport using AxisArrays, Unitful
const axes = Base.axes  # due to possible conflict with AxisArrays.axes

include("utils.jl")
include("core.jl")
include("binning.jl")

export frogtrace, freqdomain, frequencymarginal, delaymarginal
export sample, integrate, zeropad

end # module
