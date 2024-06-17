# IIR2FIR.jl
using DSP, Plots, FFTW

function sampled_freq_FIR(numTaps::Int, thisFilter::FilterCoefficients, evalFreqs::AbstractArray, dboffset::Float64=1.0)
    HofZ = freqz(thisFilter, evalFreqs, 2*evalFreqs[end])
    sampled_freq_FIR(numTaps, HofZ, evalFreqs, dboffset)
end

function sampled_freq_FIR(numTaps::Int, HofZ::AbstractArray, evalFreqs::AbstractArray, dboffset::Float64=1.0)
    numSamp = length(evalFreqs)
    HofZReal = real.(HofZ)
    HofZImag = imag.(HofZ)

    # bandpass behavior
    # phase = 0 at DC
    numSamp = length(evalFreqs)
    radpersamp = -π/2 * (2*numSamp - 1) / (numSamp)

    # Set the phase
    for j in 1:numSamp
        arg = radpersamp * j
        HofZImag[j] = HofZReal[j] * sin(arg)
        HofZReal[j] = HofZReal[j] * cos(arg)
    end

    # make symmetric in freq
    HofZReal = vcat(HofZReal, reverse(HofZReal[2:end]))
    HofZImag = vcat(HofZImag, reverse(-1 .*HofZImag[2:end]))

    # The Fourier Transform requires the center freq bins to be 0 for BPF
    HofZReal[numSamp] = 0.0
    HofZImag[numSamp] = 0.0


    R = ifft(HofZReal + im.*HofZImag)

    FIRCoef = real.(R[numSamp:numSamp+numTaps-1]) # zero out imaginary part
    return db2amp(dboffset) * FIRCoef
end

if false
    Fs = 44100
    f = range(0, Fs/2, length=4192)
    h = digitalfilter(Bandpass(90, 5000, fs=Fs), Butterworth(4))

    hi = impz(h, 4192)
    hback = PolynomialRatio(hi, [1])
    hsamp = PolynomialRatio(sampled_freq_FIR(1024, h, f), [1])
    hlevin = PolynomialRatio(lslevin_FIR(1024, freqz(h, f, Fs), f), [1])

    plot(f, pow2db.(abs.(freqz(h, f, 44100))), xscale=:log10, xlim=(10, 20000), ylims=(-50, 10), label="Orig Filter")
    plot!(f, pow2db.(abs.(freqz(hback, f, 44100))), xscale=:log10, xlim=(10, 20000), label="impz")
    plot!(f, pow2db.(abs.(freqz(hsamp, f, 44100))), xscale=:log10, xlim=(10, 20000), label="oversamp")

    plot(f, phasez(h, f), xscale=:log10, xlim=(10, 20000), label="Orig Filter")
    plot!(f, phasez(hback, f), xscale=:log10, xlim=(10, 20000), label="impz")
    plot!(f, phasez(hsamp, f), xscale=:log10, xlim=(10, 20000), label="oversamp")
end

function lslevin_FIR(numTaps::Int, HofZ::AbstractArray, evalFreqs::AbstractArray,
                    weights::AbstractArray=fill(1.0, length(evalFreqs)), dboffset::Float64=1.0)
    # complex Least Squares FIR Filter Design using Levinson's algorithm

    # evalFreqs - frequency grid (0 <= Fs/2)
    Om = (π / evalFreqs[end]) .* evalFreqs

    HofZReal = real(HofZ)
    HofZImag = imag(HofZ)
    a = fill(0.0, numTaps)
    b = fill(0.0, numTaps)

    # Setup vectors for quadratic objective
    dvec = copy(HofZ)
    evec = fill(1.0 + 0.0*im, length(seqFreq))
    e1 = exp.(im * Om)
    for i in 1:numTaps
        a[i] = transpose(weights) * real(evec)
        b[i] = transpose(weights) * real(dvec)
        evec .*= e1
        dvec .*= e1
    end

    a ./= length(numTaps)
    b ./= length(numTaps)

    return levin(a, b)
end

function levin(a, b)
    # solves system of complex linear equations toeplitz(a) * x = b using Levinson
    # a - first row of positive definite Hermetian Toeplitz matrix
    # b - right hand side vector
    n = length(a)
    t = 1
    alph = a[1]
    x = b[1]/alph

    for i in 1:n-1
        k = -(dot(a[i+1:-1:2],t))/alph
        t = vcat(t, 0) + k*reverse([conj(t); 0], dims=1)
        alph *= (1 - abs2(k))
        k = (b[i+1] - dot(a[i+1:-1:2], x)) / alph
        x = vcat(x, 0) + k*reverse(conj(t), dims=1)
    end
    return x
end
