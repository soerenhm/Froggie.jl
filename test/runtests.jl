using Froggie, Test
using AxisArrays, Unitful

@testset "Froggie" begin

  # @testset "Binning" begin
  #   include("binning.jl")
  # end

  @testset "Utils" begin
    include("utils.jl")
  end

end
