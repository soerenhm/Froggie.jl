has_delay_axis(A::AxisArray) = :t in axisnames(A)
has_frequency_axis(A::AxisArray) = :ω in axisnames(A)
has_wavelength_axis(A::AxisArray) = :λ in axisnames(A)

is_frog_trace(A::AxisArray) = has_delay_axis(A) && (has_wavelength_axis(A) || has_frequency_axis(A))


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
    # Perhaps I shouldn't call `vec2range` here...
    t = vec2range(enforce_unit(DefaultUnits.time, delays)), 
    λ = vec2range(enforce_unit(DefaultUnits.wavelength, wavelengths))
  )
end

"""
    freqdomain(data::AxisArray) -> AxisArray

Maps the wavelength axis to an angular frequency axis.
"""
function freqdomain(data::AxisArray)
  # check that data has λ axis
  has_wavelength_axis(data) || throw(ArgumentError("`data` does not containt a wavelength axis."))
  
  λ = data[Axis{:λ}] |> collect
  ω = wavelen2angfreq.(λ)

  # To convert the spectral content from the spatial domain to the time domain,
  # multiply the trace  by multiply by dλ/dω ∝ 1/ω²
    
  ξ = (ω[1] ./ ω).^2 |> ustrip
  new_data = reverse(data .* ξ', dims=2)
  ω = reverse(ω)

  AxisArray(
    new_data;
    t = AxisArrays.axes(data, Axis{:t}),
    ω = enforce_unit(DefaultUnits.frequency, ω),
  )
end


function frequencymarginal(A::AxisArray)
  nothing
end

function delaymarginal(A::AxisArray)
  nothing
end