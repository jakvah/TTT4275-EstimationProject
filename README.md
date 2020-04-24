# TTT4275-EstimationProject
Code to run simulations for TTT4275 Estimation, Detection and Classification project @ NTNU

### Maximum Likelihood Estimator (MLE) and CRLB
The code in `main.py` does the following once (!):

**a)**

Estimates the frequency of a complex sinusoidal signal embedded in complex gaussian white noise. It uses the maximum argument of the corresponding FFT as an MLE estimator for the frequency of the signal. The MLE estimator for the phase if found by multiplying the frequnecy component of the complex sinusiodal with the Fourier transform evaluated at the estimated frequency. The global variables to be changed depending on the simulation requirements are:

- Signal to noise ratio (SNR) `SNR_db` measured in decibel
- The FFT length `k` which is given as 2^k

**b)**

As long FFT estimates are undesirable in practice, an FFT of length 1024 data points is fine tuned by using a numerical search method for the frequency that minimizes the mean square error (MSE) between its own FFT and the aformentioned computed FFT.

### Automation
In order to perform and analyze larger simulations `auto.py` has been provided. It simulates point **a)** from above, but iterates through FFTs of lengths 2^{10,12,14,16,18,20}. For each FFT length iteration the signal is simulated with SNRs of {-10,0,10,20,30,40,50,60} [dB]. The amount of iterations are 1000. The results are saved to an Excel Spreedsheet file. 