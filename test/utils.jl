@testset "Padding" begin
  A = [1 2; 3 4]
  B = [0 0; 1 2; 3 4; 0 0]
  C = [0 1 2 0; 0 3 4 0]
  D = [0 0 0 0; 0 1 2 0; 0 3 4 0; 0 0 0 0]
  @test zeropad(A, 1, 1) == B
  @test zeropad(A, 1, 2) == C
  @test zeropad(A, 1) == D
  E = AxisArray(A, x = 1:2, y = [-1,2])
  @test zeropad(E, 1, Axis{:x}) == AxisArray(B, x=0:3, y=[-1,2])
  @test zeropad(E, 1, Axis{:y}) == AxisArray(C, x=1:2, y=[-4,-1,2,5])
  @test zeropad(E, 1) == AxisArray(D, x=0:3, y=[-4,-1,2,5])
end

@testset "Integration" begin

end
