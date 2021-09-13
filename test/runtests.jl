using Froggie, Test
using AxisArrays, Unitful

@testset "Froggie" begin

  @testset "Binning" begin
    include("binning.jl")
  end

end
