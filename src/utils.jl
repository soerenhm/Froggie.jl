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
    zeropad(a, N [dim::Int]) -> Array
    zeropad(a::AxisArray, N [dim::Int]) -> AxisArray

Pads the beginning and end `a` with `N` zeroes. Padding is performed
only along dimension `dim` if it is specified.
"""
zeropad(a, N, dim::Int) = copyto!(_zeropad_dim(a, N, dim)..., a, CartesianIndices(a))
zeropad(a, N) = copyto!(_zeropad_all(a, N)..., a, CartesianIndices(a))
zeropad(a::AxisArray, N, dim::Int) = AxisArray(zeropad(a.data, N, dim), expand_axes(a, N, dim)...)
zeropad(a::AxisArray, N, ax) = zeropad(a, N, axisdim(a, ax))
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

expand_axis(ax::Axis{name,T}, N) where {name,T} = Axis{name}(_expand_range(ax.val, N))
expand_axes(a::AxisArray, N) = map(ax -> expand_axis(ax, N), AxisArrays.axes(a))
expand_axes(a::AxisArray, N, dim::Int) = ntuple(n -> n == dim ? expand_axis(AxisArrays.axes(a,n),N) : AxisArrays.axes(a,n), Val(ndims(a)))

function _expand_range(x, N)
  dx1 = x[2]-x[1]
  dx2 = x[end]-x[end-1]
  [range(x[1]-N*dx1, step=dx1, length=N); x; range(x[end]+dx2, step=dx2, length=N)]
end
_expand_range(x::AbstractRange, N) = (dx = step(x); x[1]-N*dx:dx:x[end]+N*dx)


"""
    sample(a, x::AbstractArray, dim::Int, N) -> Array
    sample(a::AxisArray, ax::Axis, N) -> AxisArray
    sample(a::AxisArray, dim::Int, N) -> AxisArray

Samples dimension `dim` of `a` on `N` points distributed 
uniformly over `x`.
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
sample(a::AxisArray, ax::Type{<:Axis}, N) = sample(a, axisdim(a, ax), N)

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
    integrate(y, x::AbstractVector, dim::Int) -> Array
    integrate(y::AxisArray, dim::Int) -> AxisArray
    integrate(y::AxisArray, ax) -> AxisArray

Integrates specified axis in `y` over `x`.
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