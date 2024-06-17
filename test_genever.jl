#test_genever.jl

using Plots
using LinearAlgebra
using Random
using DSP

# Requires:
include("filterbank.jl") #mkbank
include("tf_diff_mat.jl") #tf_mat, tf_diff_mat
include("geneverStructs.jl")
include("mkEff.jl")
include("essx.jl")
include("GEneVer.jl")
include("IIR2FIR.jl")

# tests and scratch work for genever.jl

fs = 44100
lenSeq = 8192
minFreq = 90
maxFreq = 5000
numSpkr = 11

# Random.seed!(1234)
using MAT
file = matopen("/Users/matthew.baron/Perforce/audio-engineering/DSP/project_apollo/P1aFA4ArrayDesign/polars/ApolloP1a_FA4_20181031_raws.mat")
x1 = read(file, "x1")
close(file)
x1 = reshape(x1, (8192, 72, 1, 1, numSpkr))
ang = collect(-180:5:175)
dIdx = findall(.&((ang .< 30), (ang .> -30)))
uIdx = findall(.|((ang .> 30), (ang .< -30)))
# desired[samples, azimuths, altitudes, microphones, speakers]
desired = x1[:,dIdx, :, :, :]
undesired = x1[:,uIdx,:,:,:]
interaural = fill(1.0, (4097, 36, 1, 2, numSpkr))

seqFreq = freq(periodogram(desired[:,1,1,1,1], fs=fs, nfft=2*(lenSeq-1)))

weightEfficiency = fill(1.0, lenSeq)
weightEffort = 1e9*rand(8192, numSpkr)

normDesired = true
normUndesired = true
flagInterAural = false
bandstopAdd = 2
fbThreshold = -100 # db

# refIA = # can be per frequency
# weightOptIA = # can be per freq
optFilterbank = [1/3, 20, 4, 0, 0] # [octaves, order, overlap, bandstopAdd, alignFirst]


inDataStruct = inData(desired, undesired, interaural, weightEfficiency, weightEffort, normDesired, normUndesired, flagInterAural, bandstopAdd, true, 1, 1, fbThreshold)

eq,gData = genever(inDataStruct, fs)
# try out different smoothing
plot(seqFreq, amp2db.(abs.(eq[:,1])), xscale=:log10, xlims=(10, 20000), label="GEneVer out")
# plot!(seqFreq, amp2db.(essx(abs.(eq[:,1]), 1/3, "sg")), xscale=:log10, xlims=(10, 20000), label="sg smoothed")
plot!(seqFreq, amp2db.(essx(abs.(eq[:,1]), 1/6, "erg")), xscale=:log10, xlims=(10, 20000), label="erg smoothed")

# Fit these things:
eq1FIR = PolynomialRatio(sampled_freq_FIR(8192, eq[:,1], seqFreq, 6.0), [1.0])
eq1FIRsg = PolynomialRatio(sampled_freq_FIR(8192, essx(abs.(eq[:,1]), 1/3, "sg"), seqFreq), [1.0])
eq1FIRerg = PolynomialRatio(sampled_freq_FIR(8192, essx(abs.(eq[:,1]), 1/3, "erg"), seqFreq), [1.0])

plot(seqFreq, unwrap(angle.(eq[:,1])), xscale=:log10, xlim=(10, 20000), label="GEneVer Filter")
plot!(seqFreq, unwrap(phasez(eq1FIR, seqFreq)), xscale=:log10, xlim=(10, 20000), label="Direct Fit")
# plot!(seqFreq, unwrap(phasez(eq1FIRsg, seqFreq)), xscale=:log10, xlim=(10, 20000), label="SG Fit")
plot!(seqFreq, unwrap(phasez(eq1FIRerg, seqFreq)), xscale=:log10, xlim=(10, 20000), label="erg Fit")


graph_freq_resp!(eq1FIR, seqFreq, fs, "Direct Fit")
graph_freq_resp!(eq1FIRsg, seqFreq, fs, "SG Fit")
graph_freq_resp!(eq1FIRerg, seqFreq, fs, "erg Fit")

eq1LFIR = PolynomialRatio(lslevin_FIR(8192, eq[:,1], seqFreq), [1])
graph_freq_resp!(eq1FIR, seqFreq, fs, "Sampled Fit")
graph_freq_resp!(eq1LFIR, seqFreq, fs, "Levinson Fit")

graph_imp_resp(eq1FIR, 8192, "Sampled Fit")
graph_imp_resp(eq1LFIR, 8192, "Levinson Fit")

plot(seqFreq, unwrap(angle.(eq[:,1])), xscale=:log10, xlim=(10, 2000), label="GEneVer Filter")
plot!(seqFreq, unwrap(phasez(eq1FIR, seqFreq)), xscale=:log10, xlim=(10, 2000), label="Direct Sampled Fit")
plot!(seqFreq, unwrap(phasez(eq1LFIR, seqFreq)), xscale=:log10, xlim=(10, 2000), label="Direct Levinson Fit")
ylims!(-720, 720)
plot!(seqFreq, unwrap(phasez(eq1LFIRerg, seqFreq)), xscale= xlim=(10, 2000), label="ERG Levinson Fit")

ergTarg = dropdims(essx((eq[:,1]), 1/6, "erg"), dims=2)
eq1LFIRerg = PolynomialRatio(lslevin_FIR(8192, ergTarg, seqFreq), [1])
graph_freq_resp!(eq1LFIRerg, seqFreq, fs, "Levin Fit erg Smooth")
graph_imp_resp(eq1LFIRerg, 8192, "Levin Fit erg")

function graph_freq_resp(b, frs, fs, label)
    r = freqz(b, frs, fs)
    plot(frs, amp2db.(abs.(r)), xaxis=:log10, xlim=(10, fs/2), label=label)
end

function graph_freq_resp!(b, frs, fs, label)
    r = freqz(b, frs, fs)
    plot!(frs, amp2db.(abs.(r)), xaxis=:log10, xlim=(10, fs/2), label=label)
end

function graph_periodogram(signal, fs, label)
    psd = periodogram(signal, fs=fs)
    plot!(freq(psd), pow2db.(power(psd)), label=label)
end

function graph_imp_resp(b, n::Int, label)
    i = impz(b, n)
    plot(i, label = label)
end

function graph_imp_resp!(b, n::Int, label)
    i = impz(b, n)
    plot!(i, label=label)
end

## FIR fit the resulting ZPGs
zpgbank = gData.zpgbank

plot(seqFreq, amp2db.(sum(map(x-> abs.(freqz(x, seqFreq, fs)), zpgbank[:,1]))), xscale=:log10, xlims=(20, 20000))
ylims!(-40, 10)
# plot!(seqFreq, amp2db.(sum(map(x-> abs.(freqz(x, seqFreq, Fs)), zpgbank[:,2]))), xscale=:log10, xlims=(20, 20000))

# ZPG to FIR one filter
prBank = map(x -> PolynomialRatio((impz(x, 8192)), [1]), zpgbank[:,1])
plot!(seqFreq, amp2db.(sum(x-> abs.(freqz(x, seqFreq, fs)), prBank)), label="sum of impz")

prBankSum = PolynomialRatio(sum(map(x -> x.b, prBank)), prBank[1].a)
plot!(seqFreq, amp2db.(abs.(freqz(prBankSum, seqFreq, fs))), label="sum of FIR")

## FIR Fit the bandpass and bandstop independently
bptargresp = sum(map(x-> freqz(x, seqFreq, Fs), zpgbank[:,1]))

bp1FIR = PolynomialRatio(sampled_freq_FIR(8192, bptargresp, seqFreq, 6.0), [1.0])
plot!(seqFreq, amp2db.(abs.(freqz(bp1FIR, seqFreq, Fs))))

bp1FIRlev = PolynomialRatio(lslevin_FIR(8192, bptargresp, seqFreq), [1.0])
plot!(seqFreq, amp2db.(abs.(freqz(bp1FIRlev, seqFreq, Fs))))


## filter the responses with the appropriate filters, and export to .mat file to visualize the polars
eqFIR = Array{Float64}(undef, 8192, numSpkr)
for i in 1:numSpkr
    eqFIR[:, i] = lslevin_FIR(8192, dropdims(essx((eq[:,1]), 1/6, "erg"), dims=2), seqFreq)
end

y1 = Array{Float64}(undef, 8192, 72, numSpkr)
for i in 1:numSpkr
    y1[:,:,i] = filt(db2amp(23) * eqFIR[:,i], [1], x1[:,:,1,1,i])
end
file = matopen("ApolloLego_002.mat", "w")
write(file, "y1", y1)
close(file)
