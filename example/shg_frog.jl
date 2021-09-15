using Revise
using Froggie, FFTW
using JSON, DelimitedFiles
import PyPlot as plt; plt.pygui(true)


datafile = joinpath(@__DIR__, "shg-frog_001.dat")
metafile = joinpath(@__DIR__, "shg-frog_001.json")

data = readdlm(datafile)
meta = open(JSON.parse, metafile)

t = Froggie.center(meta["Data"]["delays"] .|> float) * u"fs"
λ = (meta["Data"]["wavelengths"] .|> float) * u"nm"
b = meta["Data"]["background"] .|> float

trace = frogtrace(data .- b', t, λ)
trace = trace[λ=400u"nm"..480u"nm"]

freqtrace_ = freqdomain(trace)
freqtrace = zeropad(binalong(freqtrace_, Axis{:ω}, 128), 128)
freqtrace.data[findall(freqtrace .< 0)] .= 0.0


A = frequencymarginal(freqtrace)
B = FFTW.fftshift(FFTW.ifft(A))
s = sqrt.(B)

α = 0.09
β = 0.425
γ = 1.0

r = similar(s)
fill!(r, 0)
for n in length(s)>>1-1:length(s)-1
  s₊ = real(s[n+1]) >= 0 ? s[n+1] : -s[n+1]
  Δ₀ = (s₊, -s₊) .- r[n]
  Δ₁ = Δ₀ .- (r[n]-r[n-1])
  Δ₂ = Δ₁ .- (r[n]-r[n-1] - (r[n-1]-r[n-2]))
  ϵ = @. α*abs2(Δ₀) + β*abs2(Δ₁) + γ*abs2(Δ₂)
  i = argmin(ϵ)
  r[n+1] = i == 1 ? s[n+1] : -s[n+1]
end



fig = plt.figure()
plt.subplot(2,1,1)
plt.pcolormesh(trace[Axis{:λ}].|>ustrip, trace[Axis{:t}].|>ustrip, trace)
plt.colorbar()
plt.subplot(2,1,2)
plt.pcolormesh(freqtrace[Axis{:ω}].|>ustrip, freqtrace[Axis{:t}].|>ustrip, freqtrace)
plt.colorbar()