#mkEff.jl
using DSP
# Input is a 3xNxnumSpkr array, in which
# x = listE[:,1,1] represents [lowFreq, highFreq, weightInterval]
function mkEff(listE, fs, lenSeq)
    # t, numSeg, numSpkr = size(listE)
    # if t != 3
    #     error("Input array doesn't have correct first dimension shape")
    # end



    # filterSetup = remez(lenSeq, definedBands, Hz=fs, maxiter=50),
    # filterResult =

    return listE

end
