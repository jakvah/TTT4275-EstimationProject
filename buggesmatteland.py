import numpy as np

def GenerateNoise(signal, SNR_dB): # desiredSNR [dB]; signal is an array with complex values
   n = np.zeros((len(signal),1), dtype=complex)
   # SNR = var_signal / var_n
   snr = 10.0**(SNR_dB/10.0) # Desired linear SNR
   var_signal = signal.var() # Measure power of signal
   var_n = var_signal / snr # Calculate required noise power for desired SNR
   if (var_n == 0): # In case of a null signal
       var_n = snr
   e = np.random.normal(0, np.sqrt(var_n*2.0)/2.0, size=(len(signal), 2)).view(np.complex) # Generate noise with calculated power
   for i in range(0, len(signal)):
       n[i,0] = n[i,0] + e[i][0]
   print(10.0*np.log10((np.var(signal))/(np.var(n))))
   return n

