# matplotlib.use("Agg") # Use Agg if you want to draw figs on numerical servers
import matplotlib.pyplot as plt

from numpy import double, arange, floor, zeros, pi, sin, shape, real, imag, absolute, angle, exp, random, argmax, sqrt, arctan2, mean
from scipy import fftpack, signal

class TimeSeries(object):
    """
    Member variables
    ----------------
    data   : original time series data
    fs     : sampling frequency
    N      : length of data
    DeltaF : frequency resolution
    NyN    : Nyquist frequency
    time   : Real time
    freq   : Frequency sequence which is used to show the power spectol

    fft_amp   : Amplitude
    fft_phase : Phase
    fft_ps    : Power Spectol
    fft_sine  : sin component of fourier data
    fft_cosine: cosin component of fourier data

    Member functions
    ----------------
    fft : do spectol analysis

    Sample
    ------
    To know the usage of this class,
    just run the following command:
    python FourieAnalysisClass.py
    """
    def __init__(self,data,fs,realdata=True,Trend=True):
        """
        Parameters
        ----------
        data : time series data what you want to do spectol analysis
        fs   : sampling frequency
        """
        # if np.shape(data) != '(1,)':
        #     print "error : data must be 1 dimensional data!!"
        #     sys.exit(1)
        self.data = data
        self.fs   = fs
        self.N    = len(data)
        self.DeltaF = double(fs)/double(self.N)
        self.NyN    = int(floor(0.5*self.N)+1)
        self.time   = arange(self.N)/fs
        # self.freq = np.fft.fftfreq(self.NyN, self.DeltaF)
        self.freq = fftpack.fftfreq(self.NyN, self.DeltaF)
        # self.freq = arange(self.NyN)*self.DeltaF

        if Trend:
            self.FilterTrend()

        if realdata:
            self.rfft()

            self.corr = real(fftpack.ifft(self.fft_ps))
            self.tau = self.time[:self.NyN]*2.0
            self.corr /= self.corr[0]
            # self.corr *= 2.0/self.N
        else:
            self.fft()
            
            self.corr = real(fftpack,ifft(self.fft_ps))
            self.tau = self.time[:self.NyN]*2.0
            self.corr /= self.corr[0]

    def rfft(self):
        data_fft = fftpack.rfft(self.data)
        self.fft_amp = zeros(self.NyN)
        self.fft_phase=zeros(self.NyN)
        N=self.N
        self.fft_amp[0] = data_fft[0]/double(N)
        self.fft_phase[0] = 0.0
        if self.N%2 == 0:
            for i in range(1,self.NyN-1):
                self.fft_amp[i] = sqrt( data_fft[2*i-1]**2 + data_fft[2*i]**2)/double(N)
                self.fft_phase[i] = arctan2( data_fft[2*i-1], data_fft[2*i])*180.0/pi
        else:
            for i in range(1,self.NyN-1):
                self.fft_amp[i] = sqrt( data_fft[2*i-1]**2 + data_fft[2*i]**2)/double(N)
                self.fft_phase[i] = arctan2( data_fft[2*i-1], data_fft[2*i])*180.0/pi
            
        self.fft_ps = self.fft_amp*self.fft_amp
        _exp=self.fft_amp*exp(-1.0j*self.fft_phase*pi/180.0)
        self.fft_cosine = real(_exp)
        self.fft_sine   = imag(_exp)
        return self.fft_amp, self.fft_phase

    def fft(self):
        data_fft = fftpack.fft(self.data)
        self.fft_amp = zeros(self.NyN)
        self.fft_phase=zeros(self.NyN)
        N=self.N
        for i,_data_fft in enumerate(data_fft[:self.NyN]):
            if (i==0):
                self.fft_amp[i] = real(_data_fft)/double(N)
                self.fft_phase[i] = 0.0
            elif i==self.N*0.5:
                self.fft_amp[i] = absolute(_data_fft)*(1.0/double(N))
                self.fft_phase[i] = 0.0
            else:
                self.fft_amp[i] = absolute(_data_fft)*(2.0/double(N))
                self.fft_phase[i] = angle(_data_fft)*180.0/pi

        self.fft_ps = self.fft_amp*self.fft_amp
        _exp=self.fft_amp*exp(-1.0j*self.fft_phase*pi/180.0)
        self.fft_cosine = real(_exp)
        self.fft_sine   = imag(_exp)
        return self.fft_amp, self.fft_phase

    def ShowPowerSpectorum(self):
        plt.figure()
    
        plt.subplot(211)
        plt.plot(self.time, self.data)
        plt.xlabel("time(sec)")
        plt.ylabel("data")
        
        plt.subplot(212)
        plt.loglog(self.freq, self.fft_ps)
        plt.xlabel("frequency(Hz)")
        plt.ylabel("Power Spectol")
        
        plt.show()

        # ft_fft = lf.LeastSquareOptimize(self.freq, self.fft_ps, 'power')
        # ft_fft.optimize_mask_data(0.0001, -1.0, 2, 100)
        # [a,b] = ft_fft.GetOptParams()
        # ft_fft.ShowFittingFunc('log')

    def SavePowerSpectorum(self,figname):
        plt.figure()
    
        # plt.subplot(211)
        # plt.plot(self.time, self.data)
        # plt.xlabel("time(sec)")
        # plt.ylabel("data")
        
        # plt.subplot(212)
        plt.loglog(self.freq, self.fft_ps)
        plt.xlabel("frequency(Hz)")
        plt.ylabel("Power Spectol")
        
        plt.savefig(figname)
        plt.close()

        # ft_fft = lf.LeastSquareOptimize(self.freq, self.fft_ps, 'power')
        # ft_fft.optimize_mask_data(0.0001, -1.0, 2, 100)
        # [a,b] = ft_fft.GetOptParams()
        # ft_fft.ShowFittingFunc('log')

    def CalcAutoCorr(self, corr):
        # data_fft = np.fft.fft(self.data)
        
        corr = fftpack.ifft(self.fft_ps)
        # self.fft_amp = zeros(self.NyN)
        # self.fft_phase=zeros(self.NyN)
        # N=self.N
        # for i,_data_fft in enumerate(data_fft[:self.NyN]):
        #     if (i==0):
        #         self.fft_amp[i] = real(_data_fft)/double(N)
        #         self.fft_phase[i] = 0.0
        #     elif i==self.N*0.5:
        #         self.fft_amp[i] = absolute(_data_fft)*(1.0/double(N))
        #         self.fft_phase[i] = 0.0
        #     else:
        #         self.fft_amp[i] = absolute(_data_fft)*(2.0/double(N))
        #         self.fft_phase[i] = angle(_data_fft)*180.0/pi

        # self.fft_ps = self.fft_amp*self.fft_amp
        # _exp=self.fft_amp*exp(-1.0j*self.fft_phase*pi/180.0)
        # self.fft_cosine = real(_exp)
        # self.fft_sine   = imag(_exp)
        # return self.fft_amp, self.fft_phase

    def FilterTrend(self):
        ave = mean(self.data, dtype=np.float)
        self.data -= ave

if __name__== '__main__':

    print('Sample Usage')
    T  = 10000 # define the time steps of data
    fs = 100.0   # define the sampling frequency

    data = random.random(T)
    ts = TimeSeries(data, fs)

    plt.figure()

    plt.subplot(211)
    plt.plot(ts.time, ts.data)
    plt.xlabel("time(sec)")
    plt.ylabel("data")
    plt.title('Random Noise')

    plt.subplot(212)
    plt.loglog(ts.freq, ts.fft_ps)
    plt.xlabel("frequency(Hz)")
    plt.ylabel("Power Spectol")

    plt.show()

    freq = 0.001 # frequency of sine wave
    data = sin( (2.0*pi)*arange(0,T,0.01)*freq ) + 0.5*sin( (2.0*pi)*arange(0,T,0.01)*freq*2.0 )
    ts = TimeSeries(data, fs)
    indx = argmax(ts.fft_ps)
    print(ts.freq[indx], ts.fft_ps[indx])
    ts.ShowPowerSpectorum()

    # plt.figure()
    
    # plt.subplot(211)
    # plt.plot(ts.time, ts.data)
    # plt.xlabel("time(sec)")
    # plt.ylabel("data")
    # plt.title('Sin wave with frequency '+str(freq))

    # plt.subplot(212)
    # plt.loglog(ts.freq, ts.fft_ps)
    # plt.xlabel("frequency(Hz)")
    # plt.ylabel("Power Spectol")

    # plt.show()
