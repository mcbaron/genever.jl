#filterbank.jl
using DSP, LinearAlgebra
# Implement mkbank
# [filterbank, bpFreqs, bandStop] = mkbank(Fl, Fh, Fs, lenFreq, octaves, order, overlap)
# Create a filterbank with zero phase linkwitz-riley
# fl - low freq
# fh - high freq (if fl=fh or fh is empty only one filter is returned around fl)
# Fs - sample rate (default 44100)
# lenFreq - frequency length, assumes to Nyquist (4097)
# octaves - octaves to smooth
# order - order of butterworth (doubled as LR)
# overlap - number of times to overlap filters

# filterbank - output
# bpFreqs - array of low & high freqs for each filter
# bandStop - lowpass & high pass outside filterbank

function mkbank(Fl, Fh; Fs=44100, lenFreq=4097, octaves=1/3, order=5, overlap=4)
    filtMethod = Butterworth(order)

    Fn = Fs/2
    freqRes = Fn/(lenFreq - 1)
    freqList = collect(0:freqRes:Fs/2)

    if octaves == 0 # no smoothing, linear case
        Fl = floor(Fl/freqRes)*freqRes
        Fh = ceil(Fh/freqRes)*freqRes
        FlInd = findfirst(Fl, freqList)
        FhInd = findfirst(Fh, freqList)
        numTaps = FhInd - FlInd + 1
        bpFreqs = [ freqList[FlInd:FhInd]; freqList[FlInd:FhInd] ]'
        bandStop = [digitalfilter(Lowpass(Fl, fs=Fs), filtMethod), digitalfilter(Highpass(Fh, fs=Fs), filtMethod)]
        bandStop = bandStop .* bandStop
        filterbank = [fill(0.0,(numTaps, FlInd-1)), Matrix{Float}(I,numTaps, numTaps), fill(0.0, (numTaps, lenFreq-FhInd))]
        return filterbank, bpFreqs, bandStop
    end

    # Center Freq list
    Fh = Fl*2^(ceil(log2(Fh/Fl)/octaves)*octaves)
    freqArray = Fl*2 .^ (0.0:octaves/overlap:log2(Fh/Fl))
    if isempty(Fh) || Fl == Fh # This is a center frequency
        freqMin = Fl*2^(-octaves / 2)
        freqMax = Fl*2^(octaves / 2)
        freqArray = [freqMin freqMax]
    end

    filterbank = Array{ZeroPoleGain}(undef, length(freqArray) - overlap)
    bpFreqs = Array{Float64}(undef, length(freqArray) - overlap, 2)
    for hpInd = 1:length(freqArray) - overlap
        lpInd = hpInd + overlap
        freqHp = freqArray[hpInd]
        freqLp = freqArray[lpInd]
        bp = digitalfilter(Bandpass(freqHp, freqLp, fs=Fs), filtMethod)
        bp = bp * bp
        filterbank[hpInd] = bp
        bpFreqs[hpInd,:] = [freqHp, freqLp]
    end

    #normalize filterbank - necessary?
    freqMin = freqArray[1]
    freqMax = freqArray[end]
    freqMid = 10^((log10(freqMax*freqMin))/2)
    filterbankSum = sum(map(x -> (abs.(freqz(x, freqList, Fs))), filterbank))
    normVal = filterbankSum[bisect_left(freqList, freqMid)]
    filterbank *= 1/normVal


    # bandstop calculation
    bandStop = [digitalfilter(Lowpass(freqMin, fs=Fs), filtMethod), digitalfilter(Highpass(freqMax, fs=Fs), filtMethod)]
    bandStop = bandStop .* bandStop

    return filterbank, bpFreqs, bandStop
end

# See test_filterbank for scratch space and helpful commands

function bisect_left(a, x, lo=1, hi=nothing)
    lo < 1 && throw(BoundsError(a, lo))
     hi == nothing && (hi = length(a))

     while lo < hi
         mid = (lo + hi) ÷ 2
         a[mid] < x ? lo = mid + 1 : hi = mid
     end
     return lo
 end
