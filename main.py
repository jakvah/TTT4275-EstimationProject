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
w = np.random.randn(N, 2).view(np.complex128)

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
    
    plt.figure(1)
    plt.plot(w)
    plt.title("White complex Gaussian noise")
    plt.savefig("noise.png")
    
    plt.figure(3)
    plt.plot(s)
    plt.title("Onyl signal")
    plt.savefig("signal.png")

    plt.figure(2)
    plt.title("Total sigal")
    plt.xlabel("n")
    plt.ylabel("x[n] = ")
    plt.plot(x)
    plt.savefig("total.png")    




main()

