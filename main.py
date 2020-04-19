import numpy as np
import matplotlib.pyplot as plt
import buggesmatteland as bml

# ---------- Specifications ---------- #
A = 1
SNR = 10 # In dB
SIGMA_SQUARED = (A**2)/(2*SNR)
T = 10**(-6) 
N = 513
n_0 = -256
f_0 = 10**5
omega_0 = 2*np.pi*f_0
theta = np.pi/8

k = 10
fft_length = 2**k

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
    print("The CRLB for the Omega estimator is: ", CRLB_OMEGA)
    print("The CRLB for the Theta estimator is: ",CRLB_THETA)
    
    n = np.arange(N)
    
    # White noise
    plt.figure(1)
    plt.plot(w)
    plt.title("White complex Gaussian noise")
    plt.savefig("noise.png")
    
    # Only signal
    plt.figure(2)
    plt.stem(n,s)
    plt.title("Plot of DT signal: $s[n] = e^{j(\omega n T -1)}$")
    plt.savefig("signal.png")

    # Signal + white noise
    n = np.arange(N)
    plt.figure(3)
    plt.title("Plot of DT signal: $x[n] = e^{j(\omega n T -1)} + w[n]$")
    plt.xlabel("n")
    plt.ylabel("x[n]")
    plt.stem(n,x)
    plt.savefig("stem.png")

    # Fourier transform
    FT_s = np.fft.fft(s,n = fft_length)
    FT_x = np.fft.fft(x,n = fft_length)


    print(len(x))
    print(len(s))
    print()
    print("FT_s length is ", len(FT_s))
    print("FT_x length is ", len(FT_x))
    
    
    plt.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.4)
    plt.subplot(211)
    plt.title("FFT of only signal")
    plt.stem(np.arange(len(FT_s)),FT_s)
    plt.subplot(212)
    plt.title("FFT of signal embedded in noise")
    plt.stem(np.arange(len(FT_x)),FT_x)
    plt.savefig("fft.png")    
    

    # Checking that most dominant is as is to be expected
    f = bml.findDominantFrequency(FT_s,T,fft_length)
    print("Most dominant frequency in s[in] is: ", f, " Hz")

    # Finding most dominant in total signal
    f = bml.findDominantFrequency(FT_x,T,fft_length)
    print("Most dominant frequency in x[n] is: ",f, " Hz")

main()

