#! /usr/bin/env python

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import multiprocessing as mp
import manageevent as me
import time, optspex, sys

# SPectrum Analysis, Reduction, and Calibration (SPARC) pipeline
# Written by Kevin B. Stevenson in June of 2012
import w1reduce as w1
# import w2reduce as w2
import w3lc as w3
import w4ld as w4
# import w5lc as w5

eventname = 'hat38v1'
eventdir   = '2017-05-17-w2-256'

def interactive():
    
    # w1 or w2 are the choses for "Reuciton Packages"
    # w1 is default; w2 is in production
    # eventname: user defined identifier for this project/visit
    # isplots: Flag for "if plot" and "what plots"
    #   options: 1,3,5,6,7 (bigger numbers include smaller numbers; i.e. 5 includes 1 & 3, etc)
    #       1: Flat field and background images
    #       3: Frame dependent info; i.e. 1D spectra, 2D frames w/ & w/o background subtract + extraction box
    #       5: Stuffidy stuff (for testing one frame)
    #       6: More Stuffidy stuff (for testing one frame)
    #       7: The ultimate Stuffidy stuff (for testing one frame)
    w1.reduceWFC3(eventname, isplots=3)   #Reduction + extraction
    #ev=w2.reduceWFC3(eventname, isplots=1)
        
    # White Light Curve -- because nchan = 1
    # w3 requires w1 output
    w3.lcWFC3(eventname, eventdir, 1, isplots=1, expand=1)
    
    # w5 requires w2 output
    # w5.lcWFC3   (eventname, eventdir, 1, isplots=1, expand=1)
    
    # Spectroscopic Light Curves -- because nchan > 1; nchan == 10 here -- so 10 channels
    w3.lcWFC3    (eventname, eventdir, 10, isplots=1, expand=1)
    
    for i in range(len(wave_low)):
        print(i,len(wave_low))
        # Light curve
        w3.lcWFC3   (eventname, eventdir, wave_low[i], wave_hi[i], i, isplots=1, expand=1)
        # Limb darkening
        evname = eventname + '_' + str(int(wave_low[i]*1e3)) + '_' + str(int(wave_hi[i]*1e3))
        w4.ld_driver(evname, eventdir, wave_low[i], wave_hi[i], n_param=2, isplots=True, stellarmodel='kurucz')

    # Load data for testing purposes
    ev           = me.loadevent(eventdir + '/d-' + eventname + '-w1')
    #ev           = me.loadevent(eventdir + '/d-' + evname + '-w3')
    #ev           = me.loadevent(eventdir + '/d-wa043bhp1_1125_1650-w3-WHITE')
    aux          = me.loadevent(eventdir + '/d-' + eventname + '-data')
    ev.spectra   = aux.spectra
    ev.specerr   = aux.specerr
    ev.specbg    = aux.specbg
    #ev.data      = handle['data']
    ev.data_mhdr = aux.data_mhdr
    ev.data_hdr = aux.data_hdr
    ev.mask      = aux.mask

    # Plot image
    plt.figure(2)
    plt.clf()
    plt.imshow(ev.mask[1][0][ev.ywindow[0][0]:ev.ywindow[0][1], ev.xwindow[0][0]:ev.xwindow[0][1]], origin='lower')
    plt.imshow(ev.data[32][1][ev.window[0,1]:ev.window[1,1]], origin='lower')

    # Plot spectra
    plt.figure(1)
    plt.clf()
    plt.errorbar(ev.wave[0][0], ev.spectra[32,0,2], ev.specerr[32,0,2])
    plt.vlines(wave_low,0,np.max(ev.spectra),'r','dashed')
    plt.vlines(wave_hi,0,np.max(ev.spectra),'r','dashed')

    plt.figure(3)
    plt.clf()
    plt.plot(ev.wave[0], ev.specbg[32,0])

#MAIN IS CALLED WHEN EXECUTING DIRECTLY FROM BASH
def main(args):
    length = len(args)
    #Run p6model
    if args[1] == 'p1':
        p1.reduceGMOS(*args[2:])
    #Run p7anal
    elif args[1] == 'p2':
        p2.wavecal(*args[2:])
    #Run p8tables
    elif args[1] == 'p3':
        p3.lcGMOS(*args[2:])
    #Run p9figs
    elif args[1] == 'p4':
        p4.ld_driver(*args[2:])
    else:
        print("Unrecognized function.")
    return 0

#CALLS main IN BASH MODE THEN EXITS CLEANLY
if __name__ == '__main__':
    sys.exit(main(sys.argv))

