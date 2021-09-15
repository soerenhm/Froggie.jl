center(x) = x .- (first(x)+last(x))/2

vec2range(x::AbstractVector) = range(first(x), last(x), length=length(x))

enforce_unit(unit::Unitful.FreeUnits, x::T) where T <: Number = x * unit
enforce_unit(unit::Unitful.FreeUnits, x::Q) where Q <: Quantity = uconvert(unit, x)
enforce_unit(unit::Unitful.FreeUnits, x::AbstractVector{T}) where T = enforce_unit.(Ref(unit), x)

wavelen2angfreq(λ) = uconvert(DefaultUnits.frequency, 2π*Unitful.c0 / enforce_unit(DefaultUnits.wavelength, λ))

module DefaultUnits
  using Unitful
  const time = u"fs"
  const wavelength = u"nm"
  const frequency = u"PHz"
end


"""
    sample(a, x::AbstractArray, dim::Int, N) -> Array
    sample(a::AxisArray, ax::Axis, N) -> AxisArray
    sample(a::AxisArray, dim::Int, N) -> AxisArray

Samples `a` on `N` points uniformly distributed over `x` along 
a specified axis or dimension of `a`.
"""
function sample(a, x::AbstractArray, dim::Int, N)
  @boundscheck size(a, dim) == length(x)
  Ipre = CartesianIndices(size(a)[1:dim-1])
  Ipost = CartesianIndices(size(a)[dim+1:end])
  T = promote_type(eltype(a), eltype(ustrip(x)))
  b = zeros(T, size(Ipre)..., N, size(Ipost)...)
  _sample!(b, a, x, range(first(x), last(x), length=N), Ipre, Ipost)
end

function sample(a::AxisArray, ax::Axis{name,T}, N) where {name, T}
  dim = axisdim(a, ax)
  b = sample(a.data, ax.val, dim, N)
  x = range(extrema(ax.val)..., length=N)
  axs = ntuple(n -> n == dim ? Axis{name}(x) : AxisArrays.axes(a,n), Val(ndims(a)))
  AxisArray(b, axs...)
end
sample(a::AxisArray, dim::Int, N) = sample(a, AxisArrays.axes(a, dim), N)

@noinline function _sample!(b, a, x, y, Ipre, Ipost)
  j = 2
  for i in 2:length(y)-1
    while y[i] > x[j]
      j += 1
    end
    t = (y[i]-x[j-1])/(x[j]-x[j-1])
    b[Ipre, i, Ipost] = t*a[Ipre, j, Ipost] + (1-t)*a[Ipre, j-1, Ipost]
  end
  b[Ipre, 1, Ipost] = a[Ipre, 1, Ipost]
  b[Ipre, length(y), Ipost] = a[Ipre, length(x), Ipost]
  b
end


"""
    zeropad(a, N [dim::Int]) -> Array
    zeropad(a::AxisArray, N [dim::Int]) -> AxisArray

Pads the beginning and end `a` with `N` zeroes. If `dim` is specified,
the padding is performed only along this dimension.
"""
zeropad(a, N, dim::Int) = copyto!(_zeropad_dim(a, N, dim)..., a, CartesianIndices(a))
zeropad(a, N) = copyto!(_zeropad_all(a, N)..., a, CartesianIndices(a))
zeropad(a::AxisArray, N, dim::Int) = AxisArray(zeropad(a.data, N, dim), expand_axes(a, N, dim)...)
zeropad(a::AxisArray, N) = AxisArray(zeropad(a.data, N), expand_axes(a, N)...)

function _zeropad_all(a, N)
  adims = Val(ndims(a))
  newsize = ntuple(n -> size(a,n)+2N, adims)
  Rdst = CartesianIndices(ntuple(n -> (N+1:N+size(a,n)), adims))
  zeros(eltype(a), newsize), Rdst
end

function _zeropad_dim(a, N, dim)
  adims = Val(ndims(a))
  newsize = ntuple(n -> n == dim ? size(a,n)+2N : size(a,n), adims)
  Rdst = CartesianIndices(ntuple(n -> n == dim ? (N+1:N+size(a,n)) : axes(a,n), adims))
  zeros(eltype(a), newsize), Rdst
end

function expand_axis(ax::Axis{name,T}, N) where {name,T}
  dx1 = ax.val[2]-ax.val[1]
  dx2 = ax.val[end]-ax.val[end-1]
  Axis{name}([range(ax.val[1]-N*dx1, step=dx1, length=N); ax.val; range(ax.val[end]+dx2, step=dx2, length=N)])
end

expand_axes(a::AxisArray, N) = map(ax -> expand_axis(ax, N), AxisArrays.axes(a))
expand_axes(a::AxisArray, N, dim::Int) = ntuple(n -> n == dim ? expand_axis(AxisArrays.axes(a,n),N) : AxisArrays.axes(a,n), Val(ndims(a)))