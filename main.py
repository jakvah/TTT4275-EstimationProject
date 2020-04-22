import numpy as np
import matplotlib.pyplot as plt
import buggesmatteland as bml
import math
import statistics as st

# ---------- Specifications ---------- #
A = 1.0
SNR_db = -10 # In dB

SNR_linear = 10.0**(SNR_db/10)
SIGMA_SQUARED = (A**2)/(2*SNR_linear)

T = 10**(-6) 
N = 513
n_0 = -256
f_0 = 10**5
omega_0 = 2*np.pi*f_0
theta = np.pi/8

k = 22
fft_length = 2**k

ITERATIONS = 10

# ---------- CRLB Helpers ---------- #
P = (N*(N-1)) / 2
Q = (N*(N-1)*(2*N-1)) / 6

# ---------- CRLB ---------- #

CRLB_OMEGA = (12*(SIGMA_SQUARED)) / ((A**2)*(T**2)*N*((N**2)-1))
CRLB_THETA = 12*(SIGMA_SQUARED)*((n_0**2)*N + 2*n_0*P + Q) / ((A**2)*(N**2)*((N**2)-1))

def iterate():
    # ---------- Signals ---------- #

    # White complex Gaussian noise
    wReal = np.random.normal(0, np.sqrt(SIGMA_SQUARED), size=N)
    wImag = np.random.normal(0, np.sqrt(SIGMA_SQUARED), size=N)*1j

    w = []
    for i in range(N):
        w.append(wReal[i] + wImag[i])

    # Exponential signal
    s = []
    for n in range(N):
        s.append(A*np.exp(np.complex(0,1)*(omega_0)*n*T + theta))

    # Total signal
    x = []
    for i in range(N):
        x.append(s[i] + w[i])

    # Fourier transform
    FT_s = np.fft.fft(s,n = fft_length)
    FT_x = np.fft.fft(x,fft_length)

    FT_S = np.absolute(FT_s)
    FT_x = np.absolute(FT_x)
       
    # Checking that most dominant is as is to be expected
    f_1 = bml.findDominantFrequency(FT_s,T,fft_length)
    #("Most dominant frequency in s[in] is: ", f_1/1000, " kHz")

    # Finding most dominant in total signal
    f_2 = bml.findDominantFrequency(FT_x,T,fft_length)

    return f_2

def main():    
    print("Running ",ITERATIONS, "iterations with:")
    print("SNR [dB]:",SNR_db)
    print("FFT length:",fft_length,"(2^" + str(k) + ")")
    print("Frequency:",f_0/1000,"kHz")  
    print("The CRLB for the Omega estimator is:", ((CRLB_OMEGA)/(2*np.pi))/1000,"kHz")
    print("The CRLB for the Theta estimator is:",CRLB_THETA)
    print()
    
    print(" *--------------- RESULTS ---------------*")


    error=[]
    freqs = []
    for i in range(ITERATIONS):
        f = iterate()
        err = f_0 - f
        print("Iteration nr",(i+1),": ",f/1000,"kHz. Error:",(err/1000),"kHz")
        freqs.append(f)
        error.append(err) 

    errmean=st.mean(error)
    print("Mean freq error is: ", errmean/1000, "kHz")

    errvar=st.variance(error, errmean)
    print("The variance of the freq error is: ",errvar/1000, "kHz")
     
main()