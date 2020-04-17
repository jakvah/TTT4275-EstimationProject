import numpy as np
import matplotlib.pyplot as plt
import buggesmatteland as bml

# ---------- Specifications ---------- #
A = 1
SNR = 10 # In dB
SIGMA_SQUARED = (A**2)/(2*SNR)
T = 10**(-6) 
N = 50
n_0 = -256
f_0 = 10**5
omega_0 = 2*np.pi*f_0
theta = np.pi/8

# ---------- CRLB Helpers ---------- #
P = (N*(N-1)) / 2
Q = (N*(N-1)*(2*N-1)) / 6

# ---------- CRLB ---------- #

CRLB_OMEGA = (12*(SIGMA_SQUARED)) / ((A**2)*(T**2)*N*((N**2)-1))
CRLB_THETA = 12*(SIGMA_SQUARED)*((n_0**2)*N + 2*n_0*P + Q) / ((A**2)*(N**2)*((N**2)-1))

# ---------- Signals ---------- #

# White complex Gaussian noise
wReal = np.random.normal(0, 1, size=N)
wImag = np.random.normal(0, 1, size=N)*1j

w = []
for i in range(N):
    w.append([wReal[i] + wImag[i]])

# Exponential signal
s = []
for n in range(N):
    s.append(A*np.exp(np.complex(0,1)*(omega_0)*n*T - 1))

# Total signal
x = []
for i in range(N):
    x.append(s[i] + w[i])


def main():
    n = np.arange(N)
    print("The CRLB for the Omega estimator is: ", CRLB_OMEGA)
    print("The CRLB for the Theta estimator is: ",CRLB_THETA)
    
    plt.figure(1)
    plt.plot(w)
    plt.title("White complex Gaussian noise")
    plt.savefig("noise.png")
    
    plt.figure(3)
    plt.stem(n,s)
    plt.title("Plot of DT signal: $s[n] = e^{j(\omega n T -1)}$")
    plt.savefig("signal.png")

    plt.figure(2)
    plt.title("Total sigal")
    plt.xlabel("n")
    plt.ylabel("x[n] = ")
    plt.plot(x)
    plt.savefig("total.png") 

    n = np.arange(N)
    plt.figure(4)
    plt.title("Plot of DT signal: $x[n] = e^{j(\omega n T -1)} + w[n]$")
    plt.xlabel("n")
    plt.ylabel("x[n]")
    plt.stem(n,x)
    plt.savefig("stem.png")

    # Fourier transform


    FT_s = np.fft.fft(s)
    FT_x = np.fft.fft(x)

    plt.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.4)
    plt.subplot(211)
    plt.title("FFT of only signal")
    plt.stem(n,FT_s)
    plt.subplot(212)
    plt.title("FFT of signal embedded in noise")
    plt.stem(n,FT_x)
    plt.savefig("fft.png")    

main()

