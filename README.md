Test to validate error in FFTW3.3.10 computation of the Discrete Fourier Transform.

The test computes various 1D Real to Complex DFT's of 2^n such that the values can
all be compared. If the FFTW3 values are valid then they should match the respective
output of the other 3 test DFT computations.

The resultant output is the computation of a DFT with 64, 128, 256, 512 and 1024 points
and the respective comparison made available via a CSV that can loaded as a
spreadsheet.

The respective CSV files included with the project come from running this test on
both a MacOS, x64 Windows 11 and Ubuntu builds of FFTW v3.3.10. Note that both the
64 and 128 point DFT computations are all good. But for 256, 512 and 1024 the values
computed by FFTW3 are wrongly computed when compared against the other DFT calculations.

Now I have computed the DFT using the slow direct method both on the Real to Complex and 
Complex to Complex. Also I compute the DFT a third way using a Radix-2 FFT algorithm. I
computed it via three different methods as I thought I was making mistake in my DFT
calculations NOT expecting to find the problem lying with FFTW3.

I built FFTW3.3.10 from source using the CMake project.

Stuart MacLean
Big Bang Kilonova Hypothesis www.bigbangkilonova.org
Black Omega www.blackomega.co.uk
