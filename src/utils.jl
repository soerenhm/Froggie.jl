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