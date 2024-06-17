#test_filterbank.jl

# Tests and scratch space for filterbank.jl

using Plots
Fs = 44100
Fl = 1174
Fh = 1479
filtMethod = Butterworth(20)
responseType = Highpass(1200, fs=Fs)
mylp = digitalfilter(Lowpass(Fh, fs = Fs), filtMethod)
mylp = mylp * mylp
myhp = digitalfilter(Highpass(Fl, fs = Fs), filtMethod)
myhp = myhp * myhp
f = range(10, stop=Fs/2, length=4097)
myLR = myhp * mylp

plot(f, amp2db.(abs.(freqz(myLR, f, Fs))), xscale=:log10, xlims=(10, 20000))
plot!(f, amp2db.(abs.(freqz(.5*myLR, f, Fs))))
ylims!((-40, 5))


plot(f, phasez(myLR, f./Fs), xaxis=log)
plot!(f, phasez(myButter, f./Fs), xaxis=log)

filterbank, _, __ = mkbank(Fl, Fh)

# Plot the filterbank as individual filters
plot(f, map(x -> amp2db.(abs.(freqz(x, f, Fs))), filterbank), xscale=:log10, xlims=(10, 20000))
# Plot the reconstruction of the filterbank
plot!(f, amp2db.(sum(map(x -> abs.(freqz(x, f, Fs)), filterbank))), xscale=:log10, xlims=(10, 20000))

# filterbankSum = prod(filterbank)
# plot(f, amp2db.(abs.(freqz(filterbankSum, f, Fs))), xscale=:log10, xlims=(10, 20000))

## testing init
Fl = 90
Fh = 5000
Fs = 44100
lenFreq = 4097
octaves = 1/3
order = 10
overlap = 1





# remez?
f = collect(range(10, stop=Fs/2, length=4097))
a = (collect(f) .- .5, collect(f) .+ .5)
b = reshape(collect(Iterators.flatten(a)), 20001, 2)
c = [(b[i, :]..., ) for i=1:size(b, 1)]
m = Dict(zip(c, abs.(freqz(filterbank[5], f, Fs))))
h = remez(4097, collect(f), abs.(freqz(fb[5])), Hz=Fs)
