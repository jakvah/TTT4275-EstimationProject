import numpy as np
import matplotlib.pyplot as plt
import buggesmatteland as bml
import cmath
import statistics as st
import scipy.optimize

# ---------- Specifications ---------- #
A = 1.0
SNR_db = 60 # In dB

SNR_linear = 10.0**(SNR_db/10)
SIGMA_SQUARED = (A**2)/(2*SNR_linear)

T = 10**(-6) 
N = 513
n_0 = -256
f_0 = 10**5
omega_0 = 2*np.pi*f_0
theta = np.pi/8

k = 10
fft_length = 2**k

ITERATIONS = 10

# ---------- CRLB Helpers ---------- #
P = (N*(N-1)) / 2
Q = (N*(N-1)*(2*N-1)) / 6

# ---------- CRLB ---------- #

CRLB_OMEGA = (12*(SIGMA_SQUARED)) / ((A**2)*(T**2)*N*((N**2)-1)) # In Radians^2
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
        s.append(A*np.exp(np.complex(0,1)*((omega_0)*(n + n_0)*T + theta)))
    

    # Total signal
    x = []
    for i in range(N):
        x.append(s[i] + w[i])

    # Fourier transform
    FT_x = np.fft.fft(x,fft_length)
    
    # Finding most dominant in total signal
    f_2,i = bml.findDominantFrequency(np.absolute(FT_x),T,fft_length)

    t = np.angle((np.exp(-(np.complex(0,1)*2*np.pi*f_2*n_0*T)))*FT_x[i])

    return f_2,t

def functionToBeMinimized(f_variable,initalEstimateFFT):
    # Exponential signal
    s = []
    for n in range(N):
        s.append(A*np.exp(np.complex(0,1)*((f_0)*(n + n_0)*T + theta)))

    fftGuess = np.fft.fft(s,2**10)
    
  
    return bml.meanSquareError(fftGuess,initalEstimateFFT)
    

def main():

    print("Running ",ITERATIONS, "iterations with:")
    print("SNR [dB]:",SNR_db)
    print("FFT length:",fft_length,"(2^" + str(k) + ")")
    print("Frequency:",f_0/1000,"kHz")  
    print("The CRLB for the Omega estimator is:", CRLB_OMEGA/(4*np.pi**2),"Hz^2")
    print("The CRLB for the Theta estimator is:",CRLB_THETA)
    print()
    
    print(" *--------------- RESULTS ---------------*")

    error_theta=[]
    error_f=[]
    freqs = []
    thetas = []
    for i in range(ITERATIONS):
        f,t = iterate()
        err_f = f_0 - f
        err_theta=theta-t
        freqs.append(f) # In hertz
        error_f.append(err_f) 
        thetas.append(t)
        error_theta.append(err_theta)


    print("Mean freq:",st.mean(freqs))
    errmean=st.mean(error_f)
    print("Mean freq error is: ", errmean/1000, "kHz")
    thetamean = st.mean(thetas)

    errvar_f=st.variance(error_f, errmean)
    print("The variance of the freq error is: ",errvar_f, "Hz^2")

    mean_error=st.mean(error_theta)

    print("Mean theta is:", thetamean)
    print("Mean theta error is: ", st.mean(error_theta))
    print("The variance of the phase error is: ", st.variance(error_theta, mean_error))

    bml.circlePlot(int(bml.makeAnglePositive(thetamean*180/np.pi)))
    
    print("Doing part b)")
    print()
     
    # Generate a signal with some sweet sweet noise
    # White complex Gaussian noise
    wReal = np.random.normal(0, np.sqrt(SIGMA_SQUARED), size=N)
    wImag = np.random.normal(0, np.sqrt(SIGMA_SQUARED), size=N)*1j

    w = []
    for i in range(N):
        w.append(wReal[i] + wImag[i])

    # Exponential signal
    s = []
    for n in range(N):
        s.append(A*np.exp(np.complex(0,1)*((omega_0)*(n + n_0)*T + theta)))
    

    print("S:")
    print(s)


    # Total signal
    x = []
    for i in range(N):
        x.append(s[i] + w[i])
        
    print("X:")
    print(x)

    # Fourier transform
    FT = np.fft.fft(x,2**10)
    print("lenght of initial:",len(FT))
    initalGuess = f_0 - 30000

    result = scipy.optimize.minimize(functionToBeMinimized,initalGuess,FT,method = "Nelder-Mead")
    print(result.x)



main()