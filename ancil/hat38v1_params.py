import numpy as np

# Define event object
class event_init(object):
    def __init__(self):
        
        # Unique event name
        self.eventname  = 'wasp79b'               # Name of the event; defined by user
        self.eventdir   = '2017-06-21-w1'   # Directory name; defined by user; date is defined by user
        
        # Planet parameters
        self.period     = 3.6623817         # Planetary Period
        self.perioderr  = 5.1e-06          # Planetary Period Error
        self.ephtime    = 2455545.23479     # Ephemeris timing
        
        # Variables
        self.expand     = 1                 # Expansion factor for interpolation: used in np.zoom(expand)   # Only used as 1 or 2
                                            # Useful for background subtraction
        self.diffthresh = 10                # Sigma theshold for bad pixel identification in the differential non-destructive reads (NDRs)
        self.sigthresh  = [5,5]             # Sigma threshold for bad pixel identification in the temporal access -- useful for W2 (not here)
        self.p3thresh   = 5                 # Sigma threshold for bad pixel identification in the step 3 of the Horne paper (outliers in the background subtraction)
        self.p5thresh   = 10                # Sigma threshold for bad pixel identification in the step 5 of the Horne paper (outliers in the dispersion access)
        self.p7thresh   = 10                # Sigma threshold for bad pixel identification in the step 7 of the Horne paper (outliers in the cross-dispersion access)
        self.fittype    = 'smooth'          # Generating profile along the dispersion access; 'smooth': fitting to smoothed spectrum
        self.window_len = 11                # Length scale for smoothing kernel when generating the weights for the weighted integrations (optimal spectral extraction)
        self.deg        = None              # Used for polynomial instead of smoothign function to determine the wieghts in the optimal spectral extraction
        self.bgdeg      = 0                 # Order of polynomial for background subtraction on a column by column basis (per NDR, per column)
        self.flatsigma  = 30                # Sigma threshold for flat field
        # self.filter     = 'G102'            # Label of the filter profile (code now reads the header)
        
        # Analysis inputs
        self.ncpu       =  4                # 1: for testing and figures; 8: for awesome
        self.n_spec     =  1                # Number of spectra on detector (hopefully just 1)
        self.n_stars    =  1                # Number of stars   on detector (hopefully just 1)
        # self.sep_back   = 30                # ?
        self.spec_width = 17                # Half-width of spectral extraction region (in the cross-dispersion)
        self.fitbghw    = 17                # Distance from the spectrum (in the cross-dispersion) to start background subtraction
        self.iref       = [1,2]             # Used for W2 (1 = forward scan, 2 = reverse scan)
        self.inormspec  = [-6,None]         # Used for W2 (uses last 6 spectra to normalize in a certain plot)

        # Specific gain levels for WFC3
        self.gain       = 1.                # Brightness unit is e-/sec (~early 2016) or e- (after that)
        # self.readNoise  = 20.               # data_mhdr['READNSEA'] -- data_mhdr['READNSED'] # extracted from header

        # Extraction windows
        self.xwindow    = [[100,260]]       # X-limits for Sub-array region to extract spectrum
        self.ywindow    = [[130,380]]       # Y-limits for Sub-array region to extract spectrum (used to compute background)
        # self.flatoffset = [[374, 374]]    # [[379, 379]]
        
        # Manually mask regions [specnum, colstart, colend, rowstart, rowend]
        #   Example: self.manmask = [[0,135,141,128,256]]
        # self.manmask = []
        
        # File names that contain lists of flat, bias, wavelength, and object files
        self.filebase  = None               # Basname for spectral image files
        self.flatfile  = "/Users/kristin/Documents/STScI/WASP79b-2017-03-03/cal/sedFFcube-both.fits"  # WFC3 IR Flat Field
        self.objfile   = 'all'              # Look at all files in the director; options: 'all' or 'filename.txt'
        self.objstart  = 0                  # Not needed if objfile == 'all'; otherwise determines which file to start accessing
        self.objend    = 0                  # Not needed if objfile == 'all'; otherwise determines which file to end   accessing
        self.directfile    = 'directImg.txt'# Text file with the filename for the direct image
        self.centroidguess = [280,355]      # Centroid in [y,x]
        self.centroidtrim  = 5              # Sub-array size (probably half-width) fro centroiding
        
        # File locations for science and calibration data
        self.loc_sci  = '/Users/kristin/Documents/STScI/WASP79b-2017-03-03/sci/'                                  # Science working directory with spectral images (ONLY!)
        self.loc_cal  = '/Users/kristin/Documents/STScI/WASP79b-2017-03-03/cal/'            # Direct image fits file
        self.loc_ld   = '/Users/kristin/Documents/STScI/WASP79b-2017-03-03/cal/'    # List of limb darkening
        self.horizonsfile = '/Users/kristin/Documents/STScI/WASP79b-2017-03-03/cal/2017a_hst.vec'              # convertions from UTC to TDB
        
        # File locations and names for limb darkening
        self.file_med = self.loc_ld + 'hatp38_med.txt'  # High temperature LDCs # Not working
        self.file_low = self.loc_ld + 'hatp38_low.txt'  # Mid  temperature LDCs # Not working
        self.file_hi  = self.loc_ld + 'hatp38_hi.txt'   # Low  temperature LDCs # Not working
        
        self.havecalaor = False                         # Always False

# Define auxiliary object
class aux_init(object):
    def __init__(self):
        pass
