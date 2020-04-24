import numpy as np
import matplotlib.pyplot as plt
import sys

def GenerateNoise(signal, SNR_dB): # desiredSNR [dB]; signal is an array with complex values
   n = np.zeros((len(signal),1), dtype=complex)
   # SNR = var_signal / var_n
   snr = 10.0**(SNR_dB/10.0) # Desired linear SNR
   var_signal = signal.var() # Measure power of signal
   var_n = var_signal / snr # Calculate required noise power for desired SNR
   if (var_n == 0): # In case of a null signal
       var_n = snr
   e = np.random.normal(0, np.sqrt(var_n*2.0)/2.0, size=(len(signal), 2)).view(np.complex) # Generate noise with calculated power
   for i in range(0, len(signal)):
       n[i,0] = n[i,0] + e[i][0]
   print(10.0*np.log10((np.var(signal))/(np.var(n))))
   return n

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



def circlePlot(angle):
    N = 360
    bottom = 0
    max_height = 4

    theta = np.linspace(0.0, 2 * np.pi, N, endpoint=False)
    radii = max_height*np.random.rand(N)
    width = ((2*np.pi) / N) + 0.01

    for i in range(N):
        radii[i] = 0
        if i == int(round(angle)):
            radii[i] = (max_height*1)

    ax = plt.subplot(111, polar=True)
    bars = ax.bar(theta, radii, width=width, bottom=bottom)
    ax.set_yticklabels([])
    ax.set_title("Plot av vinkelen: " + str(int(round(angle))) + " i grader:")

    # Spicy farger
    for r, bar in zip(radii, bars): 
        bar.set_facecolor(plt.cm.jet(r / 10.))
        bar.set_alpha(0.8)
    plt.savefig("angleplot.png")
    plt.show()

def makeAnglePositive(ang):
    if ang < 0:
        ang = 360 + ang
    return ang