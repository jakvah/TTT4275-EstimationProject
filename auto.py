import numpy as np
import matplotlib.pyplot as plt
import buggesmatteland as bml
import math
import statistics as st
import xlsxwriter
import sys

# ---------- These are supposed to change ---------- #
dBs = [-10,0,10,20,30,40,50,60]
lengthPowers = [10,12,14,16,18,20]

ITERATIONS = 1000

# ---------- Constants ---------- #
A = 1.0
T = 10**(-6) 
N = 513
n_0 = -256
f_0 = 10**5
omega_0 = 2*np.pi*f_0
theta = np.pi/8

FILENAME = "innafor.xlsx"


# ---------- CRLB Helpers ---------- #
P = (N*(N-1)) / 2
Q = (N*(N-1)*(2*N-1)) / 6

# ---------- CRLB ---------- #

def iterate(SIGMA_SQUARED,fft_length):
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
        s.append(A*np.exp(np.complex(0,1)*((omega_0)*n*T + theta)))

    # Total signal
    x = []
    for i in range(N):
        x.append(s[i] + w[i])

    # Fourier transform
    FT_x = np.fft.fft(x,fft_length)
    FT_x = np.absolute(FT_x)
       
    # Finding most dominant in total signal
    f_2 = bml.findDominantFrequency(FT_x,T,fft_length)

    return f_2

def main():
    wb = xlsxwriter.Workbook(FILENAME)
    ws = wb.add_worksheet()

    # Colum names
    ws.write(0, 0, "FFT Length")
    ws.write(0, 1, "SNR [dB]")
    ws.write(0, 2, "Mean estimated")
    ws.write(0, 3, "Mean error")
    ws.write(0, 4, "Variance")
    ws.write(0, 5, "CRLB [Hz^2]")


    lengthIterationIndex = 0
    for p in lengthPowers:
        fft_length = 2**p
        dataIterationIndex = 0
        
        ws.write(1 + lengthIterationIndex*len(dBs), 0, str(fft_length))
        print("Computing FFTs of length", fft_length)
        if p == 18:
                print("Now also showing SNR progress:")

        for SNR_db in dBs:
            ws.write(1 + dataIterationIndex +lengthIterationIndex*len(dBs), 1, SNR_db)
            
            
            SIGMA_SQUARED = bml.sigmaSquaredFromdB(SNR_db,A)
            CRLB_OMEGA = (12*(SIGMA_SQUARED)) / ((A**2)*(T**2)*N*((N**2)-1))
            CRLB_THETA = 12*(SIGMA_SQUARED)*((n_0**2)*N + 2*n_0*P + Q) / ((A**2)*(N**2)*((N**2)-1))

            error=[]
            freq = []
            for i in range(ITERATIONS):
                f = iterate(SIGMA_SQUARED,fft_length)

                err = f_0 - f
                freq.append(f)
                error.append(err) 

            freqmean = st.mean(freq)
            errmean=st.mean(error)
            errvar=st.variance(error, errmean)

            ws.write(1 + dataIterationIndex +lengthIterationIndex*len(dBs), 2, freqmean)
            ws.write(1 + dataIterationIndex +lengthIterationIndex*len(dBs), 3, errmean)
            ws.write(1 + dataIterationIndex +lengthIterationIndex*len(dBs), 4, errvar)
            ws.write(1 + dataIterationIndex +lengthIterationIndex*len(dBs), 5, (CRLB_OMEGA/(4*np.pi**2)))

            dataIterationIndex += 1
            if p > 16:
                print("Done with",SNR_db)
        lengthIterationIndex += 1

    wb.close()      
    print("Added iterationsdata to",FILENAME)
main()