"""
    binalong(a, nbins) -> Vector
    binalong(a, dim::Int, nbins) -> Array
    binalong(a, x::AbstractVector, nbins) -> Vector
    binalong(a, x::AbstractVector, dim::Int, nbins) -> Array
    binalong(a::AxisArray, dim::Int, nbins) -> AxisArray
    binalong(a::AxisArray, ax::Axis, nbins) -> AxisArray

Arrange the elements of `a` along dimension `dim` into `n` bins.

If coordinates vector `x` is supplied, the bins are arranged 
uniformly over `x` (instead of uniformly over the indices of `a`).

Endpoints are handled simply by assuming that `a` has the same value
outside its bounds; this means that the binned value at the endpoint 
always has the same value as the endpoint of `a`.
"""
function binalong(a, x::AbstractVector, dim::Int, nbins)
  @boundscheck size(a, dim) == length(x)
  Ipre = CartesianIndices(size(a)[1:dim-1])
  Ipost = CartesianIndices(size(a)[dim+1:end])
  T = promote_type(float(eltype(a)), float(eltype(ustrip(x))))
  b = zeros(T, size(Ipre)..., nbins, size(Ipost)...)
  _binarray!(b, a, x, range(x[1], x[end], length=nbins), Ipre, nbins, Ipost)
end
binalong(a, dim::Integer, nbins) = binalong(a, axes(a, dim), dim, nbins)
binalong(a, nbins) = binalong(vec(a), 1, nbins)
binalong(a, x::AbstractVector, nbins) = binalong(vec(a), x, 1, nbins)

function binalong(a::AxisArray, dim::Int, nbins)
  x = AxisArrays.axes(a, dim)
  b = binalong(a.data, x.val, dim, nbins)
  x′ = range(extrema(x)..., length=nbins)
  axs = ntuple(n -> n == dim ? Axis{axisnames(a)[n]}(x′) : AxisArrays.axes(a, n), Val(ndims(a)))
  AxisArray(b, axs...)
end
binalong(a::AxisArray, ax::Type{<:Axis}, nbins) = binalong(a, axisdim(a, ax), nbins)
binalong(a::AxisArray, ax::Axis, nbins) = binalong(a, typeof(ax), nbins)

function _binarray!(b, a, x, y, Ipre, n, Ipost)
  Δy = step(y)
  x1 = y1 = x[1]
  j = 1
  for i in 1:n-1
    while true
      y2 = y[i]+Δy
      x2 = x[j+1]
      if y2 > x2
        b[Ipre, i, Ipost] += a[Ipre, j, Ipost] * (x2-y1)/(x2-x1)
        x1 = y1 = x2
        j += 1
      else
        b[Ipre, i, Ipost] += a[Ipre, j, Ipost] * (y2-y1)/(x2-x1)
        y1 = y2
        break
      end
    end
  end
  b[Ipre, n, Ipost] = a[Ipre, length(x), Ipost]
  b
end