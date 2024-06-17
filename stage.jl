#stage.jl

# Stage a collection of measurements and call teiger.jl

# Do we want to use SampledSignals.jl and interface with a multi-channel audio card,
# or should we read in a directory of .wav files?
# Let's start with a directory of .wav files.

# 20190613 - since we're going to start with polar data, we'll come back to this decision


using WAV
using Glob
using MAT

file = matopen("/Users/matthew.baron/Perforce/audio-engineering/DSP/project_apollo/legoBlockArrayDesign/polars/ApolloLegoCenteredTweeter_20180517_raws.mat")
x1 = read(file, "x1")
close(file)
