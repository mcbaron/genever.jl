#genever.jl

# A targeted eigenvalue solver in Julia.

function genever(inDataStruct, fs)

    # unpack inDataStruct
    # # error checking
    # if desired.fs != undesired.fs
    #     error('Desired and Undesired containters do not have matching sample rates')
    # end
    # # etc

    # desired and undesired need to have time samples as first dimension
    irDesired = inDataStruct.desired
    irUndesired = inDataStruct.undesired
    irInterAural = inDataStruct.interaural

    # Some constants from the input data
    lenSeq = size(irDesired, 1)
    seqFreq = freq(periodogram(irDesired[:,1,1,1,1], fs=fs, nfft=2*(lenSeq-1)))
    tfFreqBW = fs/lenSeq

    # weightEfficiency extraction assign for each freq
    # weightEfficiency[frequency]
    weightEfficiency = inDataStruct.efficiency

    # weightEffort extraction assign for each freq
    # weightEffort[frequency, speaker]
    weightEffort = mkEff(inDataStruct.effort, fs, lenSeq)

    fbThreshold = db2amp(inDataStruct.fbThreshold)

    # create filterbank (linear)
    fbFilterbank, fbBPFreqs, fbBandStop = mkbank(minFreq, maxFreq, Fs=fs, lenFreq=length(seqFreq))
    fbArray = transpose(hcat(map(x -> freqz(x, seqFreq, fs), fbFilterbank)...))
    fbNumbanks = size(fbArray, 1)

    # create efficiency by filterbank
    weightEfficiencyAllFreq = fbArray .* repeat(transpose(weightEfficiency), fbNumbanks, 1)
    weightEfficiencyFb = (sum(weightEfficiencyAllFreq.^2, dims=2)./ sum(fbArray.^2, dims=2)).^ 0.5
    # this might have to transition to repeat()

    tfFreqBW = fs/lenSeq
    # create matricies and solve for X
    Xfreq = fill(0.0 + 0.0im, lenSeq, numSpkr) # X[w, spkr]
    Xzpg = Array{ZeroPoleGain}(undef, fbNumbanks, numSpkr)
    XfbVect = fill(0.0 + 0.0im, fbNumbanks, numSpkr) #Xfb[filterbank, skpr] eigenvector for debug
    bpFreqLast = [0.0, 0.0]
    for fbIdx in 1:fbNumbanks
        # get current filterbank
        fbCurr = fbArray[fbIdx, :]
        bpFreqCurr = fbBPFreqs[fbIdx,:]
        println("filterbank $(fbIdx) of $(fbNumbanks) ($(bpFreqCurr[1]) Hz - $(bpFreqCurr[2]) Hz)")

        # find freq indicies above the cutoff fbThreshold
        freqIdx = abs.(fbCurr) .>= (maximum(abs.(fbCurr)) * fbThreshold)

        # create Sd
        # frequency domain
        tfFbDesired = ir2tfmat(irDesired, fbCurr, freqIdx, fs) # apply filterbin to impulse response
        if normDesired
            tfFbDesired = tfFbDesired ./ norm(tfFbDesired)
        end
        Sd = db2amp(94).* tfFbDesired

        # create Su
        tfFbUndesired = ir2tfmat(irUndesired, fbCurr, freqIdx, fs)
        if normUndesired
            tfFbUndesired = tfFbUndesired ./ norm(tfFbUndesired)
        end
        Su = db2amp(94) .* tfFbUndesired

        if flagInterAural
            irFbInterAural = filt(fbFilterbank[fbIdx], irInterAural)
            tfFbInterAuralDiff = tf_diff_mat(irFbInterAural, fs, inDataStruct.rMicFactor, inDataStruct.weightIAFreq)
            tfFbInterAuralSum = tf_diff_mat(irFbInterAural, fs, -1, inDataStruct.weightIAFreq)
            Sd = vcat(Sd, tfFbInterAuralSum)
            Su = vcat(Su, tfFbInterAuralDiff)
        end

        Sdd = transpose(Sd)*Sd
        Suu = transpose(Su)*Su

        # Scale efficiency by norm
        Sdd_norm = norm(Sdd)
        Suu_norm = norm(Suu)
        weightMatrixEfficiency = weightEfficiencyFb[fbIdx, :] .* Matrix(I, numSpkr, numSpkr) #TODO: Check this line
        weightMatrixEfficiency = Suu_norm .* weightMatrixEfficiency

        # weight effort for out of band frequencies
        weightEffortFb = repeat(fbCurr, 1, numSpkr) .* weightEffort
        weightMatrixEffort = diagm(0=>dropdims((sum(weightEffortFb.^2, dims=1)./sum(repeat(fbCurr, 1, numSpkr), dims=1)).^ 0.5,dims=1))
        weightMatrixEffort = Suu_norm .* weightMatrixEffort

        # solve
        F = eigen(Sdd, Suu + weightMatrixEffort + weightMatrixEfficiency)
        currX = (F.vectors[:,end])

        # save vector for debug
        XfbVect[fbIdx, :] = F.vectors[:,end]

        # Assign eigenvector to filterbank
        currXwfreq = fbCurr * transpose(currX)
        currXwzpg = [fbFilterbank[fbIdx] * currX[i] for i in 1:numSpkr]

        # find the best phase to assign for this answer
        alignPhase = true
        if (alignPhase && sum(bpFreqLast) != 0.0)
            lpfreq = bpFreqLast[2]
            fLPidx = findall(x -> x >= (lpfreq-tfFreqBW/2) && x <= (lpfreq+tfFreqBW/2), seqFreq[:])
            hpfreq = bpFreqCurr[1]
            fHPidx = findall(x -> x >= (hpfreq-tfFreqBW/2) && x <= (hpfreq+tfFreqBW/2), seqFreq[:])
            # find the average phase difference over the overlap
            deltaPhase = exp(im*angle(sum(Xfreq[[fHPidx; fLPidx],:] .* conj(currXwfreq[[fHPidx; fLPidx],:]))))
            # # use phase on the zpg TODO: This needs testing
            # deltaPhase = exp(im*sum(sum(phasez(XzpgLast, [fHPidx:fLPidx]./fs) .- phasez(currXwzpg, [fHPidx:fLPidx]./fs) )))
            currXwfreq *= deltaPhase
            currXwzpg *= deltaPhase
            XfbVect[fbIdx, :] *= deltaPhase
        end
        # force first bin to closest 1/-1
        if (inDataStruct.alignFirstBand && sum(bpFreqLast) == 0.0)
            currIdxFreqBP = findall(x -> x > bpFreqCurr[1] && x < bpFreqCurr[2], seqFreq[:])
            Xwreal = real(currXwfreq[currIdxFreqBP,:])
            Xwreal = (Xwreal ./ Xwreal) .* abs.(currXwfreq[currIdxFreqBP,:])
            deltaPhase = exp(im*angle(sum(Xwreal.*conj(currXwfreq[currIdxFreqBP,:]))))
            currXwfreq *= deltaPhase
            currXwzpg *= deltaPhase
        end

        # Accumulate fb to answer
        Xfreq = Xfreq + currXwfreq
        Xzpg[fbIdx, :] = currXwzpg

        # update for next iteration
        bpFreqLast = copy(bpFreqCurr)
        XzpgLast = copy(currXwzpg)
    end

    # add in the bandstop
    if bandstopAdd == 1 # match edge values (band the same as edges)
        fbBandStopFreq = map(x -> freqz(x, seqFreq, fs), fbBandStop)
        Xfreq = Xfreq + transpose(XfbVect[1,:] * transpose(fbBandStopFreq[1]) + XfbVect[end,:] * transpose(fbBandStopFreq[2]))
    elseif bandstopAdd == 2 # force to be real
        fbBandStopFreq = map(x -> freqz(x, seqFreq, fs), fbBandStop)
        bsFirstReal = (real(XfbVect[1,:])/abs.(real(XfbVect[1,:]))) * abs.(XfbVect[1,:])
        bsLastReal = (real(XfbVect[end,:])/abs.(real(XfbVect[end,:]))) * abs.(XfbVect[end,:])
        Xfreq = Xfreq + transpose(bsFirstReal * transpose(fbBandStopFreq[1]) + bsLastReal * transpose(fbBandStopFreq[2]))
    end

    # scaling
    avg = sum(abs.(Xfreq))./ prod(size(Xfreq))
    Xfreq *= 1/avg

    # Reshape EQ matrix if necessary
    eq = Xfreq

    # Debug data
    gData = geneverData(XfbVect, fbFilterbank, Xzpg, eq)

    return eq, gData
end
