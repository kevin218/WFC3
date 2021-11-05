from pylab import *;ion()

def plot_fft(xdata, ydata, fStart = 0, normed = True):
		'''
		    >>> nData     = 1e4
		    >>> tSpacing  = 1.0 / 800.0
		    >>> xdata     = np.linspace(0.0, nData*tSpacing, nData)
		    >>> ydata     = np.sin(50.0 * 2.0*np.pi*xdata) + 0.5*np.sin(80.0 * 2.0*np.pi*xdata)
		    >>> epa.auxiliary.plot_fft(x,y,fStart = 0, normed=True)

		    >>> frequency    = np.linspace(0.0, 1.0/(2.0*tSpacing), nData/2)
		    >>> gaussNoise   = np.random.normal(0, 1, frequency.size)
		    >>> redNoise     = np.abs(np.fft.ifft(1. / frequency))
		    >>> fullNoise    = abs(np.fft.ifft(np.fft.fft(gaussNoise) / frequency))
		'''

		from scipy.fftpack import fft
		import matplotlib.pyplot as plt

		if normed:
				ydata   = ydata.copy()
				ydata  -= np.median(ydata)

		# Number of samplepoints
		nData   = ydata.size
		# sample spacing
		spacing     = np.median(np.diff(xdata))
		ydata_fft   = fft(ydata)
		frequency   = np.linspace(0.0, 1.0/(2.0*spacing), nData//2)

		plt.semilogy(frequency[fStart:], 2.0/nData * np.abs(ydata_fft[fStart:nData//2]))
		plt.grid()
		plt.show()

def plot_fft_with_window(xdata, ydata, fStart = 1, normed = True, **kwargs):
		'''
		    >>> # Number of samplepoints
		    >>> N = 600
		    >>> # sample spacing
		    >>> T = 1.0 / 800.0
		    >>> x = np.linspace(0.0, N*T, N)
		    >>> y = np.sin(50.0 * 2.0*np.pi*x) + 0.5*np.sin(80.0 * 2.0*np.pi*x)
		'''
		from scipy.fftpack import fft
		from scipy.signal import blackman
		import matplotlib.pyplot as plt

		nData       = ydata.size

		if fStart > nData/2:
				raise Exception("fStart must be less than nData/2 or there won't be anything to plot")

		spacing     = np.median(np.diff(xdata))
		ydata_fft   = fft(ydata)

		window      = blackman(nData)
		ydata_wfft  = fft(ydata*window)

		frequencies = np.linspace(0.0, 1.0/(2.0*spacing), nData//2)

		#plt.semilogy(frequencies[fStart:nData/2], 2.0/nData * np.abs(ydata_fft[fStart:nData/2]), '-b')
		plt.semilogy(frequencies[fStart:nData//2], 2.0/nData * np.abs(ydata_wfft[fStart:nData//2]), '-', label='FFT Window', **kwargs)
		#plt.legend(['FFT', 'FFT Window'])
		#plt.grid()
		plt.show()
import asyncio
loop = asyncio.get_event_loop()

def hello():
    loop.call_later(3, print_hello)

def print_hello():
    print('Hello!')
    loop.stop()
    
if __name__ == '__main__':
    loop.call_soon(hello)
    loop.run_forever()
def create_dirac_comb(nPts):
	dirac_comb = np.ones(nPts*2)
	dirac_comb[::2]      = 0
	dirac_comb[:npts//2] = 0
	dirac_comb[-npts//2:] = 0
	return dirac_comb

npts = int(100)
times      = np.arange(0,2*npts) - npts/2
dirac_comb = create_dirac_comb(npts)

figure(1);clf()
plot(times,dirac_comb)
figure(2);clf()
plot_fft(times,dirac_comb)
figure(3);clf()
plot_fft_with_window(times,dirac_comb)