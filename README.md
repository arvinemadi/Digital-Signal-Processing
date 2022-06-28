## Digital Signal Processing
Some examples of digital signal processing are added here.

### Radar Range and Velocity Simulation
Use of 1D FFT for range measurement and 2D FFT for range and velocity measurement of a simulated FMWC radar is shown.
CFAR technique is implemented on the output 2D FFT for reliable estimate of the target vs. local noise. CFAR adaptively changes the threshold limit for detection based on the noise surronding the points in 2D space of distance vs. velocity in the 2D FFT.
CFAR: https://en.wikipedia.org/wiki/Constant_false_alarm_rate
Good reference on system parameters for FMWC radars: https://training.ti.com/sites/default/files/docs/mmwaveSensing-FMCW-offlineviewing_4.pdf

Two versions of the file are here. One version is simpler to follow the algorithm, but runs slower. In the second version MATLAB vectorization is used for the calculations. Mainly for calculating the local noise over the training cell in the 2D FFT, convolution with a kernel of training cells is done.

### Use of Kalman filter for noise cancellation of PPG signal
