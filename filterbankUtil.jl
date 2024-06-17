# filterbankUtil.jl
using DSP
# take in the output of mkbank and return everything as FIR filter coefficients
# using both methods in IIR2FIR.jl (argument switch)


function mkFiniteBank(filterbank::Array{ZeroPoleGain}, bpFreqs, bandStop::ZeroPoleGain, seqFreq, firMode, numTaps)
    if firMode == "sampled"
        for i in 1:length(filterbank)
            filterbank[i] = PolynomialRatio(sampled_freq_FIR(4192, filterbank[i], seqFreq, 6.0), [1])
        end
        bandStop = PolynomialRatio(sampled_freq_FIR(4192, bandStop, seqFreq, 6.0), [1])
    elseif firMode == "levin"
        for i in 1:length(filterbank)
            binTemp = freqz(filterbank[i], seqFreq, 44100)
            filterbank[i] = PolynomialRatio(lslevin_FIR(4192, binTemp, seqFreq), [1])
        end
        bandStop = PolynomialRatio(lslevin_FIR(4192, bandStop, seqFreq), [1])
    else
        error("FIR transformation mode not supported")
    end

    return filterbank, bpFreqs, bandStop
end
