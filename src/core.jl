has_delay_axis(A::AxisArray) = :t in axisnames(A)
has_frequency_axis(A::AxisArray) = :ω in axisnames(A)
has_wavelength_axis(A::AxisArray) = :λ in axisnames(A)

is_frog_trace(A::AxisArray) = has_delay_axis(A) && (has_wavelength_axis(A) || has_frequency_axis(A))

function frog_trace_type(A::AxisArray)
  if !is_frog_trace(A)
    nothing #Invalid
  else
    if has_frequency_axis(A)
      nothing # some type
    elseif has_wavelength_axis(A)
      nothing # another type
    end
  end 
end


"""
    frogtrace(data::AbstractMatrix, delays, wavelengths) -> AxisArray

Constructs an `AxisArray` holding the FROG trace (`data`) with specified
delay (first dim) and wavelengths (second dim).

If units are not explicitly specified, delays and wavelengths are assumed
given in femtoseconds and nanometers. 
"""
function frogtrace(data::AbstractMatrix, delays::AbstractArray, wavelengths::AbstractArray)
  AxisArray(
    data; 
    t = enforce_unit(DefaultUnits.time, delays), 
    λ = enforce_unit(DefaultUnits.wavelength, wavelengths)
  )
end

"""
    freqdomain(data::AxisArray) -> AxisArray

Maps the wavelength axis to an angular frequency axis.
"""
function freqdomain(data::AxisArray)
  # check that data has λ axis
  has_wavelength_axis(data) || throw(ArgumentError("`data` does not containt a wavelength axis."))
  ω = wavelen2angfreq.(data[Axis{:λ}])

  # To convert the spectral content from the spatial domain to the time domain,
  # multiply the trace  by multiply by dλ/dω ∝ 1/ω²
  ξ = (ω[1] ./ ω).^2 |> ustrip
  new_data = reverse(data .* ξ', dims=2)
  ω = reverse(ω)

  AxisArray(
    new_data,
    AxisArrays.axes(data, Axis{:t}),
    Axis{:ω}(enforce_unit(DefaultUnits.frequency, ω)),
  )
end


frequencymarginal(A::AxisArray) = has_frequency_axis(A) ? integrate(A, Axis{:t}) : nothing

function delaymarginal(A::AxisArray)
  if has_frequency_axis(A) 
    integrate(A, Axis{:ω})
  elseif has_wavelength_axis(A)
    integrate(A, Axis{:λ})
  end
end