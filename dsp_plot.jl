using Plots
theme(:juno)
pyplot(size=(900, 300))
plot(randn(10), m=:o)

using DSP
PlottableFilter = Union{FilterCoefficients}
f = digitalfilter(
    Lowpass(10, fs=40),
    Butterworth(8)
)
@show typeof(f) <: PlottableFilter;

using RecipesBase
# @userplot FilterPlot
@recipe function func(pfs::AbstractVector{<:PlottableFilter}; fs=2,
                    w=linspace(0, pi, 500),
                    hz=zeros(Float32, 0),
                    domag=true, dophase=true, dodelay=false,
                    doimp=false, impn=500,
                    dostep=false, stepn=500,
                    desc=["$ix" for ix=1:length(pfs)]
    )
    # markershape --> :auto        # if markershape is unset, make it :auto
    # markercolor :=  customcolor  # force markercolor to be customcolor
    nplots = sum([domag,dophase,dodelay,doimp,dostep])
    grid      -->  true
    layout    -->  (nplots,1)
    legend    -->  true
    linewidth -->  2
    the_title = "Filter Response" * (length(pfs)>1 ? "s" : "")
    the_xlabel = "freq [Hz]   OR   time [sec]"

    @show desc
    # determine the frequencies based on w or hz
    if length(hz) > 0
        # set w according to hz
        w = hz ./ (fs/2) .* pi
    else
        hz = w ./ pi .* (fs/2)
    end

    for (pf, label) = zip(pfs, desc)
        # calculate plottable vectors that multiple plots may depend on
        if domag || dophase || dodelay
            h = freqz(pf, w)
        end
        if dophase || dodelay
            phase = unwrap(angle.(h))
        end

        # convenience state variables, per filter
        spi=0 # subplot index
        if domag
            spi+=1; this_label = (1==spi) ? label : ""
            mag = amp2db.(abs.(h))
            @series begin
                subplot :=  spi
                label   :=  label
                ylabel  :=  "mag [dB]"
                ylims   :=  ([maximum(mag)-70, maximum(mag)+10])
                if nplots == spi; xlabel := the_xlabel; end
                if      1 == spi; title  := the_title ; end
                hz, mag
            end
        end

        if dophase
            spi+=1; this_label = (1==spi) ? label : ""
            @series begin
                subplot -->  spi
                ylabel  -->  "phase [rad]"
                label   -->  this_label
                if nplots == spi; xlabel := the_xlabel; end
                if      1 == spi; title  := the_title ; end
                hz, phase
            end
        end

        if dodelay
            spi+=1; this_label = (1==spi) ? label : ""
            @series begin
                subplot -->  spi
                ylabel  --> "delay [sec]"
                label   -->  this_label
                if nplots == spi; xlabel := the_xlabel; end
                if      1 == spi; title  := the_title ; end
                hz[2:end], (-diff(phase) ./ diff(w) ./ fs)
            end
        end

        if doimp
            spi+=1; this_label = (1==spi) ? label : ""
            ir = impz(pf, impn)
            @series begin
                subplot :=  spi
                ylabel  --> "impulse"
                label   -->  this_label
                if nplots == spi; xlabel := the_xlabel; end
                if      1 == spi; title  := the_title ; end
                (0:impn-1)./fs, ir
            end
        end

        if dostep
            spi+=1; this_label = (1==spi) ? label : ""
            sr = stepz(pf, stepn)
            @series begin
                subplot -->  spi
                ylabel  --> "step"
                label   -->  this_label
                if nplots == spi; xlabel := the_xlabel; end
                if      1 == spi; title  := the_title ; end
                (0:stepn-1)./fs, sr
            end
        end
    end # for pf=pfs
end
RecipesBase.plot(pf::PlottableFilter; kwargs...) = plot([pf]; kwargs...)
