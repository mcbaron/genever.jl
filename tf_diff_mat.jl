# tf_diff_mat.jl

# create matrix of transfer functions for M mics by N speakers
# ir[samples, azimuths, altitudes, microphones, speakers]
function tf_mat(ir, fs)
    numSamp, numAz, numAt, numMic, numSpkr = size(ir)
    sfft = 2*(numSamp - 1)
    tf = mapslices(x -> power(periodogram(x, fs=fs, nfft=sfft)).^.5, ir, dims=1)

    return reshape(tf, (numSamp*numAz*numAt*numMic, numSpkr))
end


function ir2tfmat(ir, filt, freqIdx, fs)
    numSamp, numAz, numAt, numMic, numSpkr = size(ir)
    sfft = 2*(numSamp - 1)
    tf = mapslices(x -> power(periodogram(x, fs=fs, nfft=sfft)).^.5, ir, dims=1)
    # tf in linear units
    numFreq = sum(freqIdx)
    # filter the frequency of the input signal with the filter subband
    tfFilt = tf[freqIdx,:,:,:,:] .* filt[freqIdx]

    return reshape(tfFilt, (numFreq*numAz*numAt*numMic, numSpkr))
end

# Create matrix of transfer functions for N speakers by M/2 mics (l-r)
#
# ir[samples, azimuths, altitudes, microphones, speakers]
# rMicFactor - allow for non-zero ILD/ITD reference
# freqWeight - per frequency, per speaker weight
function tf_diff_mat(ir, fs, rMicFactor, freqWeight)
    numSamp, numAz, numAt, numMic, numSpkr = size(ir)
    sfft = 2*(numSamp - 1)
    tf = mapslices(x -> power(periodogram(x, fs=fs, nfft=sfft)).^.5, ir, dims=1)

    S = freqWeight .* ((tf[:, :, :, 1, :]) .- rMicFactor .* (tf[:, :, :, 2, :]))
    return reshape(S, (numSamp*numAz*numAt, numSpkr))
end
