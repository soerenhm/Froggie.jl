# I think I need to change how I implemented binning
# What we really care about is resampling to a uniform grid
# (since we're working with Fourier transforms)

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


"""
    integrate(y, x::AbstractVector, dim::Int) -> Array
    integrate(y::AxisArray, dim::Int) -> AxisArray
    integrate(y::AxisArray, ax) -> AxisArray

Integrates the dimension `dim` of `y` over `x`.
"""
function integrate(y, x::AbstractVector, dim::Int)
  @boundscheck size(y, dim) == length(x)
  T = typeof(first(y) * first(ustrip(x)))
  sz = ntuple(n -> n == dim ? 1 : size(y,n), Val(ndims(y)))
  r = Array{T}(undef, sz...)
  # Last element of `x` is handled by ignoring it
  # This is inconsitent with how I implemented binalong...
  CI = CartesianIndices(ntuple(n -> n == dim ? Base.OneTo(size(y,n)-1) : axes(y,n), Val(ndims(y))))
  integratedim!(r, y, CI, ustrip(diff(x)), dim)
end

function integrate(y::AxisArray, dim::Int)
  x = AxisArrays.axes(y, dim)
  z = integrate(y.data, x.val, dim)
  # The following is type unstable (because dimension is reduced by one I suppose)
  axs = filter(ax -> ax != x, AxisArrays.axes(y))
  AxisArray(reshape(z, length.(axs)...), axs...)
end
integrate(y::AxisArray, ax) = integrate(y, axisdim(y, ax))

@noinline function integratedim!(r, y, CI, Δx, dim)
  fill!(r, 0)
  J = last(CartesianIndices(r))
  for I in CI
    r[min(I,J)] += Δx[I[dim]] * y[I]
  end
  r
end