# geneverStructs.jl

struct filterbankGenever
    octaves
    order
    overlap
    minFreq
    maxFreq
    Fs
    lenFreq

    filterbank
    bpFreqs
    bandstop

    array
    numBanks

    bandstopAdd
    alignFirst
end


struct geneverData
    XfbVect
    filterbank
    zpgbank
    eq
end

struct inData
    desired
    undesired
    interaural

    efficiency
    effort

    # flags
    normDesired
    normUndesired
    flagInterAural
    bandstopAdd
    alignFirstBand

    rMicFactor
    weightIAFreq
    fbThreshold
end
