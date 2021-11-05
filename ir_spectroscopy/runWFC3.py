#! /usr/bin/env python

import numpy as np
import scipy as sp
import pyfits as pf
import matplotlib.pyplot as plt
import multiprocessing as mp
import manageevent as me
import time, optspex, sys

# WFC3 Data Reduction Pipeline
# Written by Kevin B. Stevenson in June of 2012
#import w1reduce as w1
import w2reduce as w2
#import w3lc as w3
import w4ld as w4
import w5lc as w5

from importlib import reload
reload(w2)

eventname = 'kt011bph2'
eventdir   = '2021-01-05-w2'

def interactive():
    #w1.reduceWFC3(eventname, isplots=1)   #Reduction + extraction
    ev=w2.reduceWFC3(eventname, eventdir, isplots=3)
    # G141
    #wave_low = [1.125]
    #wave_hi  = [1.65]
    #w3.lcWFC3   (eventname, eventdir, 1, isplots=1, expand=1)
    #ev=w5.lcWFC3   (eventname, eventdir, 1, isplots=1, expand=1)
    #ev=w5.lcWFC3   (eventname, eventdir, 15, isplots=1, expand=1)
    # G102
    # White
    ev=w5.lcWFC3   (eventname, eventdir, 1, 0,0, wmin=0.85, wmax=1.18, isplots=1, expand=1)#, offset=0)
    # Spec - 15 chan
    ev=w5.lcWFC3   (eventname, eventdir, 15, 0,0, wmin=0.85, wmax=1.18, isplots=1, expand=1)#, offset=0)
    # Spec - 1 chan
    wave_low = np.array([0.85])
    wave_hi  = np.array([1.18])
    # Spec - 15 chan
    wmin=0.85
    wmax=1.18
    nchan=15
    binsize     = (wmax - wmin)/nchan
    wave_low    = np.round([i for i in np.linspace(wmin, wmax-binsize, nchan)],3)
    wave_hi     = np.round([i for i in np.linspace(wmin+binsize, wmax, nchan)],3)
    #Spec - 28 chan
    #wave_low = np.round([i for i in np.arange(1.10995, 1.62,0.0186)],3)
    #wave_hi  = np.round([i for i in np.arange(1.12855, 1.64,0.0186)],3)
    for i in range(len(wave_low)):
        print(i,len(wave_low))
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

    return
