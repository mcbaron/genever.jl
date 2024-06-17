# GenEVer.jl

GeneralizedEigenVector beamforming

This is an older project, circa 2019 written in the most up-to-date Julia at the time, which means that it's out of date now. 

I'll come back to it and probably port it to python sooner than later. 

## General idea:
The general idea of this project is contained in generver.jl, which is the main function. 
We take in raw polar noise recordings of the device in question, indexed by angle of incidence and which speaker in the device is active. 

The idea is to define a polar pattern for the device by selecting angles at which energy is desired and undesired, posing the problem as a Rayleigh quotient and using the eigenvectors to optimize the quotient. 
This eigenvector solution is calculated for each band of a filterbank, with the idea of utilizing the filterbank reconstruciton to enable backing out the FIR filters to apply to the audio signals which are routed to each driver in the device. 

