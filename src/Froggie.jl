module Froggie


using AxisArrays, Unitful, DelimitedFiles, JSON

const axes = Base.axes  # due to possible conflict with AxisArrays.axes


greet() = print("Hello World!")

end # module
