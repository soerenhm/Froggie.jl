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