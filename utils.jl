using FFTW
using DSP

function fft2out!(out::Array{T}, s_fft::Vector{Complex{T}}, nfft::Int, r::Real, onesided::Bool, offset::Int=0) where T
    m1 = convert(T, 1/r)
    n = length(s_fft)
    if onesided
        m2 = convert(T, 2/r)
        out[offset+1] += abs(s_fft[1])*m1
        for i = 2:n-1
            @inbounds out[offset+i] += abs(s_fft[i])*m2
        end
        out[offset+n] += abs(s_fft[end])*ifelse(iseven(nfft), m1, m2)
    else
        if n == nfft
            for i = 1:length(s_fft)
                @inbounds out[offset+i] += abs(s_fft[i])*m1
            end
        else
            # Convert real FFT to two-sided
            out[offset+1] += abs(s_fft[1])*m1
            @inbounds for i = 2:length(s_fft)-1
                v = abs2(s_fft[i])*m1
                out[offset+i] += v
                out[offset+nfft-i+2] += v
            end
            out[offset+n] += abs(s_fft[n])*m1
            if isodd(nfft)
                out[offset+n+1] += abs(s_fft[n])*m1
            end
        end
    end
    out
end

# Evaluate a window function at n points, returning both the window
# (or nothing if no window) and the squared L2 norm of the window
compute_window(::Nothing, n::Int) = (nothing, n)
function compute_window(window::Function, n::Int)
    win = window(n)::Vector{Float64}
    norm2 = sum(abs2, win)
    (win, norm2)
end
function compute_window(window::AbstractVector, n::Int)
    length(window) == n || error("length of window must match input")
    (window, sum(abs2, window))
end

function localperiodogram(s::AbstractVector{T}; onesided::Bool=eltype(s)<:Real,
                     nfft::Int=nextfastfft(length(s)), fs::Real=1,
                     window::Union{Function,AbstractVector,Nothing}=nothing) where T<:Number
    onesided && T <: Complex && error("cannot compute one-sided FFT of a complex signal")
    nfft >= length(s) || error("nfft must be >= n")

    win, norm2 = compute_window(window, length(s))
    if nfft == length(s) && win == nothing && isa(s, StridedArray)
        input = s # no need to pad
    else
        input = zeros(fftintype(T), nfft)
        if win != nothing
            for i = 1:length(s)
                @inbounds input[i] = s[i]*win[i]
            end
        else
            copyto!(input, s)
        end
    end

    s_fft = T <: Real ? rfft(input) : fft(input)
    return fft2out!(zeros(fftabs2type(T), onesided ? (nfft >> 1)+1 : nfft),
                         s_fft, nfft, fs*norm2, onesided)
end
