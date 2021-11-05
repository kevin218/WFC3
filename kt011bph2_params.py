import numpy as np

# Define event object
class event_init:
    def __init__(self):

        # Unique event name
        self.eventname  = 'kt011bph2'
        #self.eventdir   = '2016-12-11-sw15-bg15-wl31-thresh10'
        self.eventdir   = '2021-01-05-w2'

        # Planet parameters
        self.period     = 4.7362083     #Colon+ (2020)
        self.perioderr  = 0.0000040
        self.ephtime    = 2457483.43047

        # Variables
        self.expand     = 1
        self.diffthresh = 10
        self.sigthresh  = [5,5]
        self.p3thresh   = 5
        self.p5thresh   = 10
        self.p7thresh   = 10
        self.fittype    = 'smooth'
        self.window_len = 31
        self.deg        = 3
        self.bgdeg      = 0
        self.flatsigma  = 50
        self.driftzoom  = [2,2]

        # Analysis inputs
        self.ncpu       =  1
        self.n_spec     =  1
        self.n_stars    =  1
        #self.sep_back   = 30
        self.spec_width = 100
        self.fitbghw    = 100
        self.iref       = [2,3]
        self.inormspec  = [-6,None]

        # Specific gain levels for WFC3
        self.gain       = 1.     #Brightness unit is already e-/sec
        #self.readNoise  = 20.    #data_mhdr['READNSEA'] -- data_mhdr['READNSED']

        # Extraction windows
        self.xwindow    = [[150,360]]
        self.ywindow    = [[ 50,500]]
        self.flatoffset = [[374, 374]]  #[[379, 379]]
        # Manually mask regions [specnum, colstart, colend, rowstart, rowend]
        #self.manmask = [[0,135,141,128,256]]

        # File names that contain lists of flat, bias, wavelength, and object files
        self.filebase  = None
        self.flatfile  = "/Users/stevekb1/Documents/data/WFC3/ancil/flats/WFC3.IR.G102.flat.2.fits"
        #self.flatfile  = "/home/kevin/Documents/UChicago/data/ancil/flats/sedFFcube-both.fits"
        self.objfile   = 'all'  #'Orbit1a.txt'
        self.objstart  = 0
        self.objend    = 0
        self.directfile    = 'directImg.txt'
        self.centroidguess = [216,244]  #[188,150]    #centroid in [y,x]
        self.centroidtrim  = 5

        # File locations for science and calibration data
        self.loc_sci  = '/Users/stevekb1/Documents/data/WFC3/KELT11b-G102/sci/'
        self.loc_cal  = '/Users/stevekb1/Documents/data/WFC3/KELT11b-G102/cal/'
        self.loc_ld   = '/Users/stevekb1/Documents/data/limbDarkening/kelt11/'
        self.horizonsfile = self.loc_cal + '2020_hst.vec'
        self.leapdir  = '/Users/stevekb1/Documents/code/POET/POET/code/ancil/leapseconds/'
        # File locations and names for limb darkening
        self.file_med = self.loc_ld + 'kelt11_med.txt'
        self.file_low = self.loc_ld + 'kelt11_low.txt'
        self.file_hi  = self.loc_ld + 'kelt11_hi.txt'

        self.havecalaor = False
        return

# Define auxiliary object
class aux_init:
    def __init__(self):
        return
