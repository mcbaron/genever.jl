# essx.jl
using FFTW
# energy smoothing applied to a transfer function
# input is assumed to be in linear units
# default to 1/3rd octave rectangular window
# return values are also linear units
function essx(tf, octaves=1/3, method="erg")
    d3d = size(tf)
    numFreq = d3d[1]
    if length(d3d) > 1
        numSpkr = d3d[2]
    else
        numSpkr = 1
    end

    if method == "complex"
        y = fill(0.0 + 0.0*im, numFreq, numSpkr)
        for s in 1:numSpkr
            z = real(ifft(tf[:,s].^2))
            shift = argmax(abs.(z))
            shift -= 1
            fz1 = fft(circshift(z, -shift))
            fz2 = localSmooth(fz1[1:Int(numFreq/2+1),:],octaves)
            fz2 = vcat(fz2, conj(fz2[Int(numFreq/2):-1:2,:]))
            z2 = circshift(real(ifft(fz2)), shift)
            y[:,s] = fft(z2)
        end

    elseif method == "sg" #Savitzky-Golay smoothing
		sg = SG(20, 3)# half-window, polynomial degree
		# TODO: derive these numbers
		y = fill(0.0, numFreq, numSpkr)
		for s in 1:numSpkr
			y[:,s] = apply_SG_filter(sg[:,1], tf[:,s].^2)
		end
    else # method is erg
        y = localSmooth(tf.^2, octaves)
    end

    y = .^(y, 0.5) # back to linear units
    return y

end


function localSmooth(tf, octaves)
    d3d = size(tf)
    numFreq = d3d[1]
    if length(d3d) > 1
        numSpkr = d3d[2]
    else
        numSpkr = 1
    end
    # integrate & prepend a zero
    fx1 = vcat(transpose(fill(0.0, numSpkr)), cumsum(tf, dims=1))

    k = 2^(octaves/2)
    i = collect(0:numFreq-1)

    # bounds
    i1 = Int.(round.(i./k) .- 1)
    i2 = Int.(clamp.(round.(i.*k), 0, numFreq-2))

    # inverse bandwidth
    d = 1 ./(i2 .- i1)
    # calculate and return
    return (fx1[2 .+ i2,:].-fx1[2 .+ i1,:]).*d
end


function vandermonde(halfWindow::Int, polyDeg::Int,T::Type=Float64)
    @assert halfWindow>=0
    @assert polyDeg>=0

    x=T[i for i in -halfWindow:halfWindow]
    n = polyDeg+1
    m = length(x)
    V = Array{T}(undef, m, n)

    for i = 1:m
        V[i,1] = T(1)
    end
    for j = 2:n
        for i = 1:m
            V[i,j] = x[i] * V[i,j-1]
        end
    end

    return V
end

function apply_SG_filter(filter::StridedVector,signal::StridedVector)
    @assert isodd(length(filter))
    halfWindow = round(Int,(length(filter)-1)/2)
    padded_signal =
	    [signal[1]*ones(halfWindow);
         signal;
         signal[end]*ones(halfWindow)]

    filter_cross_signal = conv(filter[end:-1:1], padded_signal)
    return filter_cross_signal[2*halfWindow+1:end-2*halfWindow]
end


function SG(halfWindow::Int, polyDeg::Int,T::Type=Float64)
    @assert 2*halfWindow>polyDeg

    V=vandermonde(halfWindow,polyDeg,T)
    Q,R = qr(V)
	R = vcat(R, fill(0.0, size(V,1) - size(R,1), size(R, 2)))
    sg = R\transpose(Q)

    for i in 1:size(sg,1)
        sg[i,:]*=factorial(i-1)
    end
# CAVEAT: returns the transposed matrix
    return transpose(sg)
end
