import numpy as np
import matplotlib.pyplot as plt
import sys

def meanSquareError(list1,list2):
    if len(list1) != len(list2):
        print("List lengths don't match. Exiting")
        sys.exit()

    total = 0
    for i in range(len(list1)):
        total += (list1[i] - list2[i])**2

    return total/len(list1)
    
# Returns most dominant frequency of FFT in hertz
# fft must be numpy array, absolute value
def findDominantFrequency(fft,samplingPeriod,fftLength):
    maxVal = max(fft)
    maxIndex = np.where(fft == maxVal)
    maxIndex = maxIndex[0][0]

    f = maxIndex * (1/(samplingPeriod*fftLength))
    return f, maxIndex

def sigmaSquaredFromdB(SNR_db,A):
    SNR_linear = 10.0**(SNR_db/10)
    SIGMA_SQUARED = (A**2)/(2*SNR_linear)
    return SIGMA_SQUARED

