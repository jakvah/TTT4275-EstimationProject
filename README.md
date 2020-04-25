# TTT4275-EstimationProject
Code to run simulations for TTT4275 Estimation, Detection and Classification project @ NTNU

### Maximum Likelihood Estimator (MLE) and CRLB
The code in `main.py` does the following once (!):

**a)**

Estimates the frequency of a complex sinusoidal signal embedded in complex gaussian white noise. It uses the maximum argument of the corresponding FFT as an MLE estimator for the frequency of the signal. The MLE estimator for the phase if found by multiplying the frequnecy component of the complex sinusiodal with the Fourier transform evaluated at the estimated frequency. The global variables to be changed depending on the simulation requirements are:

- Signal to noise ratio (SNR) `SNR_db` measured in decibel.
- The FFT length `k` which is given as 2^k.
- Amount of iterations `ITERATIONS` to run.

To easliy change the simualtion variables listed above you can pass them as arugments if you are running the script from the command line. 
```
$~ python3 main.py [SNR (in decibel)*] [Radix-2 power of FFT length] [Iterations]
```

E.g if you want to simulate 100 iterations with an SNR of -10 dB using FFTs of length 2^10, you run:

```
$~ python3 main.py -10 10 100
```
**b)**

As long FFT estimates are undesirable in practice, an FFT of length 1024 data points is fine tuned by using a numerical search method for the frequency that minimizes the mean square error (MSE) between its own FFT and the aformentioned computed FFT.

### Automation
In order to perform and analyze larger simulations `auto.py` has been provided. It simulates point **a)** from above, but iterates through FFTs of lengths 2^{10,12,14,16,18,20}. For each FFT length iteration the signal is simulated with SNRs of {-10,0,10,20,30,40,50,60} [dB]. The amount of iterations are 1000. The results are saved to an Excel Spreadsheet file. 

### Math functions

The necessary mathematical calculations are done in a seperate module named `buggesmatteland.py` It is imported by both `main.py` and `auto.py` and must be included in the same folder.