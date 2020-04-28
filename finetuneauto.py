import numpy as np
import matplotlib.pyplot as plt
import buggesmatteland as bml
import cmath
import statistics as st
import scipy.optimize
import sys

# ---------- Signal specifications ---------- #
A = 1.0
SNR_db = 30
SNR_linear = 10.0**(SNR_db/10)
SIGMA_SQUARED = (A**2)/(2*SNR_linear)
print("Running with SNR:", SNR_db, "dB")
T = 10**(-6) 
N = 513
n_0 = -256
f_0 = 10**5
omega_0 = 2*np.pi*f_0
theta = np.pi/8

ITERATIONS = 100
k = 10
fft_length = 2**k
# ---------- CRLB Helpers ---------- #
P = (N*(N-1)) / 2
Q = (N*(N-1)*(2*N-1)) / 6

# ---------- CRLB ---------- #

CRLB_OMEGA = (12*(SIGMA_SQUARED)) / ((A**2)*(T**2)*N*((N**2)-1)) # In Radians^2
CRLB_THETA = 12*(SIGMA_SQUARED)*((n_0**2)*N + 2*n_0*P + Q) / ((A**2)*(N**2)*((N**2)-1))

# ---------- Computes MLE of the frecuency and phase ---------- #
def computeMLE():
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
    
    # Finding most dominant frequency in total signal
    f_2,i = bml.findDominantFrequency(np.absolute(FT_x),T,fft_length)

    t = np.angle((np.exp(-(np.complex(0,1)*2*np.pi*f_2*n_0*T)))*FT_x[i])

    return f_2,t

# ---------- Function to be minimzed in part 1B) ---------- #
def functionToBeMinimized(f_variable):  
    # ---------- Generate a FFT to be used in 1B) ---------- #

    # White complex Gaussian noise
    gwReal = np.random.normal(0, np.sqrt(SIGMA_SQUARED), size=N)
    gwImag = np.random.normal(0, np.sqrt(SIGMA_SQUARED), size=N)*1j

    gw = []
    for i in range(N):
        gw.append(gwReal[i] + gwImag[i])

    # Exponential signal
    gs = []
    for n in range(N):
        gs.append(A*np.exp(np.complex(0,1)*((omega_0)*(n + n_0)*T + theta)))
        
    # Total signal
    gx = []
    for i in range(N):
        gx.append(gs[i] + gw[i])

    gFFT = np.fft.fft(gx,2**10)
    gf = bml.findDominantFrequency(np.absolute(gFFT),T,2**10)

    f_var_sliced = f_variable[0]
    
    # -------- Exponential signal without noise -------- #
    s = []
    for n in range(N):
        s.append(A*np.exp(np.complex(0,1)*((2*np.pi*f_var_sliced)*(n + n_0)*T + theta)))

    fftGuess = np.fft.fft(s,2**10)

    mse = bml.meanSquareError(np.absolute(fftGuess),np.absolute(gFFT))

    return mse
    

def main():
    finetunes = []
    errors = []
    for i in range(ITERATIONS):
        result = scipy.optimize.minimize(functionToBeMinimized,100000,method = "Nelder-Mead")
        finetunes.append(result.x[0])
        errors.append(f_0 - result.x[0])
        print(i)
    
    meanfinetunes = st.mean(finetunes)
    meanerror = st.mean(errors)
    errvar = st.variance(errors,meanerror)

    print("Average finetuned frequency:",meanfinetunes)
    print("Average finetuned error:",meanerror)
    print("Finetuned error variance:",errvar)
    


    # Plotting MSE    
    mse = []    
    t = [1,2]
    for f in range(60000,140000,100):
        t[0] = f
        mse.append(functionToBeMinimized(t))

    plt.figure(2)
    plt.title("MSE")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Mean Square Error")
    plt.plot(np.arange(60000,140000,100),mse)
    plt.savefig("Bilder/mse.png")
        

main()