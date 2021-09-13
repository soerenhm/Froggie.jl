using Revise
using Froggie, AxisArrays, Unitful, JSON, DelimitedFiles
import PyPlot as plt; plt.pygui(true)


datafile = joinpath(@__DIR__, "shg-frog_001.dat")
metafile = joinpath(@__DIR__, "shg-frog_001.json")

data = readdlm(datafile)
meta = open(JSON.parse, metafile)

t = Froggie.center(meta["Data"]["delays"] .|> float) * u"fs"
λ = (meta["Data"]["wavelengths"] .|> float) * u"nm"
b = meta["Data"]["background"]

trace = frogtrace(data .- b', t, λ)
trace = trace[λ=400u"nm"..480u"nm"]

freqtrace_ = freqdomain(trace)
freqtrace = binalong(freqtrace_, Axis{:ω}, 128)



fig = plt.figure()
plt.subplot(2,1,1)
plt.pcolormesh(trace[Axis{:λ}].|>ustrip, trace[Axis{:t}].|>ustrip, trace)
plt.colorbar()
plt.subplot(2,1,2)
plt.pcolormesh(freqtrace[Axis{:ω}].|>ustrip, freqtrace[Axis{:t}].|>ustrip, freqtrace)
plt.colorbar()