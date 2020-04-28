# TTT4275-EstimationProject
Code to run simulations for TTT4275 Estimation, Detection and Classification project @ NTNU

### Maximum Likelihood Estimator (MLE) and Carmer Rao Lower Bound (CRLB)
The code in `main.py` does the following once (!):

**a)**

Estimates the frequency of a complex sinusoidal signal embedded in complex gaussian white noise. It uses the maximum argument of the corresponding FFT as an MLE estimator for the frequency of the signal. The MLE estimator for the phase is found by multiplying the frequency component of the complex sinusiodal with the Fourier transform evaluated at the estimated frequency. The CRLB for both frequency and phase are computed for comparison with the MLE estimates. The global variables to be changed depending on the simulation requirements are:

- Signal to noise ratio (SNR) `SNR_db` measured in decibel.
- The FFT length `k` which is given as 2^k.
- Amount of iterations `ITERATIONS` to run.

To easliy change the simualtion variables listed above you can pass them as arguments if you are running the script from the command line. 
```
$~ python3 main.py [SNR (in decibel)] [Radix-2 power of FFT length] [Iterations]
```

E.g if you want to simulate 100 iterations with an SNR of -10 dB using FFTs of length 2^10, you run:

```
$~ python3 main.py -10 10 100
```

If you choose to simply execute the program without arugments or from an IDE like Pycharm where arguments cannot be passed, it will run with its default settings; SNR of -10 dB, 100 iterations and FFTs lengths of 1024. These can be changed in the script, if desired.


**b)**

As long FFT estimates are undesirable in practice, an FFT of length 1024 data points is finetuned by using a numerical search method for the frequency that minimizes the mean square error (MSE) between its the FFT of a signal with no noise and the specified frequency and the aformentioned computed FFT.

### Automation
In order to perform and analyze larger simulations `auto.py` has been provided. It simulates point **a)** from above, but iterates through FFTs of lengths 2^{10,12,14,16,18,20}. For each FFT length iteration the signal is simulated with SNRs of {-10,0,10,20,30,40,50,60} [dB]. The amount of iterations are 1000. The results are saved to an Excel spreadsheet file named "SimulationResults.xlsx" 

`finetuning.py` servers a similar purpose for point **b)**. It finetunes the estimate and calculates the mean frequency, mean error and variance over 100 iterations. These values are computed for one SNR level and are just printing to the cmd line, not save in an Excel spreadsheet 
### Math functions

The necessary mathematical calculations are done in a seperate module named `buggesmatteland.py` It is imported by both `main.py` and `auto.py` and must be included in the same folder.