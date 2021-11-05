#! /usr/bin/env python

import numpy as np
import scipy.interpolate as spi
import matplotlib.pyplot as plt
import manageevent as me
import sys, time
import hst_scan as hst
import scipy.ndimage.interpolation as spni
from astropy.io import fits

# w3 is for light curve extraction
#   It produces both white light (nchan=1) and spectroscopic (nchan>1) light curves
#   It is dependent on the output from `w1`
def lcWFC3(eventname, eventdir, nchan, wmin=1.125, wmax=1.65, expand=1, isplots=True):
    '''
    Compute photometric flux over specified range of wavelengths
    
    Parameters
    ----------
    eventname   : Unique identifier for these data
    eventdir    : Location of save file
    nchan       : Number of spectrophotometric channels
    wmin        : minimum wavelength
    wmax        : maximum wavelength
    expand      : expansion factor
    isplots     : Set True to produce plots
    
    Returns
    -------
    None
    
    History
    -------
    Written by Kevin Stevenson      June 2012
    
    '''
    
    # Load saved data
    # An event is an instance of an object
    #   loadevent: load the saved files from storage
    #       container for all of the data 
    #           event: small data structures (i.e. light curves)
    #           aux:   large data structures (i.e. image cubes, etc)
    #
    #       i.e. event.BJDTDB : returns the time array
    #            aux.spectra  : 1D spectra per NDR
    #            aux.specerr  : 1D spectra error per NDR
    #            aux.data_mhdr: master header per frame
    #            
    print("Loading saved data...")
    ev              = me.loadevent(eventdir + '/d-' + eventname + '-w1')
    aux             = me.loadevent(eventdir + '/d-' + eventname + '-data')
    ev.spectra      = aux.spectra
    specerr         = aux.specerr
    data_mhdr       = aux.data_mhdr
    
    # Determine wavelength bins
    binsize     = (wmax - wmin)/nchan   # width in bins
    wave_low    = np.round([i for i in np.linspace(wmin, wmax-binsize, nchan)],3)   # Left  edge of the wavelength bins
    wave_hi     = np.round([i for i in np.linspace(wmin+binsize, wmax, nchan)],3)   # Right edge of the wavelength bins
    # binwave     = (wave_low + wave_hi)/2. # Middle of wavelength bin
    
    # Increase resolution of spectra: uses np.zoom to oversample the image in a flux conserving interpolation
    if expand > 1:
        # note: ev.n_spec : number of spectra per frame :: hopefully just one
        print("Increasing spectra resolution...")
        # ev.spectra.shape[3] : wavelength (dispersion) dimension
        hdspectra = np.zeros((ev.n_files,ev.n_spec,ev.n_reads-1,expand*ev.spectra.shape[3]))    # hdspectra : high definition spectra
        hdspecerr = np.zeros((ev.n_files,ev.n_spec,ev.n_reads-1,expand*ev.spectra.shape[3]))    # hdspecerr : high definition spectra error
        hdwave    = np.zeros((ev.n_img,ev.n_spec,expand*ev.spectra.shape[3]))                   # hdwave    : high definition wavelength array
        
        # This is the 'zoom' step
        for n in range(ev.n_spec): # per spectrum on the image (n_spec == 1 for WFC3)
            # This operates over all 1D stellar spectrum (over time) at once
            hdspectra[:,n]  = spni.zoom(ev.spectra[:,n],zoom=[1,1,expand])
            hdspecerr[:,n]  = spni.zoom(specerr[:,n],zoom=[1,1,expand])*np.sqrt(expand)
        for m in range(ev.n_img): # n_img : number of direct images
            # Some visits have a new wavelength solution per orbit
            for n in range(ev.n_spec):
                hdwave[m,n] = spni.zoom(ev.wave[m][n],zoom=expand)
        
        # Store high defition spectra
        ev.spectra  = hdspectra
        specerr     = hdspecerr
        ev.wave     = hdwave
    
    # Correct for drift, if calculated
    if hasattr(ev, 'drift_model') and ev.drift_model != None:
        # Correct for drift :: map the motion of the spectrum across the detector
        #   provides higher precision on the wavelength solution
        print('Applying drift correction... (Old stare-mode version, may not work)')
        # Staring Mode Operations
        # ev.drift_model is defined in `w1`
        nx = ev.spectra.shape[3] # number of pixels in the wavelength direction # FINDME: CHANGED the SHAPE[2] to a SHAPE[3]
        for m in range(ev.n_files): # iterate over time
            for n in range(ev.n_spec): # iterate over number of spectra (ideally == 1)
                spline          = spi.UnivariateSpline(np.arange(nx), ev.spectra[m,n], k=3, s=0) # Compute the spline for the shift
                ev.spectra[m,n] = spline(np.arange(nx)+ev.drift_model[n,m])  # Shifts the spectrum
        # finished Stare-mode operations
    elif ev.detector == 'IR':
        # This is for Scanning Mode
        #Calculate drift over all frames and non-destructive reads
        print('Applying drift correction...')
        # hst.drift_fit calculates the drift in the 1D spectra :: Does a cross correlation
        #   ev      : class with the data
        #   preclip : left edge of spectrum
        #   preclip : right  edge of spectrum
        ev.drift, ev.drift_model, ev.goodmask = hst.drift_fit(ev, preclip=0, postclip=None, width=5*expand, deg=2, validRange=11*expand)
        
        # Correct for drift
        if ev.n_reads > 2:
            # Throw away the first NDR -- it's bad
            print('WARNING: Marking all first reads as bad.')
            istart = 1
        else:
            print('Using first reads.')
            istart = 0
        
        # Apply the Drift correction (fancy word for spline)
        nx = ev.spectra.shape[3] # number of pixels in the wavelength direction 
        for m in range(ev.n_files): # iterate over time
            for n in range(ev.n_spec): # iterate over number of spectra (ideally == 1)
                for p in range(istart,ev.n_reads-1): 
                    # Compute the spline for the shift
                    spline            = spi.UnivariateSpline(np.arange(nx), ev.spectra[m,n,p], k=3, s=0)
                    
                    # Using measured drift, not model fit
                    # `model fit` comes from the ev.drift_model
                    # `measured drift` comes from spline of order 3
                    ev.spectra[m,n,p] = spline(np.arange(nx)+ev.drift[n,m,p]) # Shifts the spectrum
        
        #Apply scan height correction
        #print('Applying scan height correction...')
        #ev.spectra  += ev.spectra[0,0]*(1-ev.scanHeight[:,:,np.newaxis,np.newaxis]/ev.scanHeight[0,0])
        #ev.spectra /= ev.scanHeight[:,:,np.newaxis,np.newaxis]
        #specerr    /= ev.scanHeight[:,:,np.newaxis,np.newaxis]
    else:
        # UVIS Stuff
        istart = 0
    
    # Assign scan direction: 0:forward vs 1:reverse
    ev.scandir = np.zeros(ev.n_files)   # Sets up all images as forward scan: modify later
    ev.n_scan0 = 0                      # Number of forward scans
    ev.n_scan1 = 0                      # Number of reverse scans
    
    try:
        scan0 = data_mhdr[0]['POSTARG2']# POSTARG2 changes for forward and reverse scan 0: first file
        scan1 = data_mhdr[1]['POSTARG2']# POSTARG2 changes for forward and reverse scan 1: first file
        for m in range(ev.n_files):
            # for every file file, check header if POSTARG2 == scan0 or scan1
            if data_mhdr[m]['POSTARG2'] == scan0:
                # Sum up number of forward scan
                ev.n_scan0 += 1
            elif data_mhdr[m]['POSTARG2'] == scan1:
                # Store scandir == 1 for reverse scanning
                ev.scandir[m] = 1
                # Sum up number of reverse scan
                ev.n_scan1 += 1
            else:
                # Something happened
                print('WARNING: Unknown scan direction for file ' + str(m) + '.')
        
        print("# of files in scan direction 0: " + str(ev.n_scan0))
        print("# of files in scan direction 1: " + str(ev.n_scan1))
    except:
        ev.n_scan0 = ev.n_files
        print("Unable to determine scan direction, assuming unidirectional.")
    
    print("Generating light curves...")
    ev.eventname2 = ev.eventname    # Store old event name
    for i in range(nchan):
        ev.wave_low  = wave_low[i]
        ev.wave_hi   = wave_hi[i]
        print("Bandpass = %.3f - %.3f" % (ev.wave_low,ev.wave_hi))
        
        # Calculate photometric flux for each spectrum
        ev.photflux     = np.zeros((ev.n_spec, ev.n_files, np.max((1,ev.n_reads-2))))   # This become the light curve (to be populated)
        ev.photfluxerr  = np.zeros((ev.n_spec, ev.n_files, np.max((1,ev.n_reads-2))))   # This become the light curve errorbars (to be populated)
        # ev.wave         = []
        for n in range(ev.n_spec): # hopefully == 1
            if ev.detector == 'IR':
                # Compute common wavelength and indices to apply over all observations
                wave = np.zeros(len(ev.wave[0][n])) # To be the globale wavelength array
                for j in range(ev.n_img): # iterate over each image
                    wave += ev.wave[j][n]
                
                wave /= ev.n_img
                
                # wave  = np.mean(ev.wave, axis=0) # FINDME: TEST LATER
                
                # index == where wave meets BOTH requirement
                # Which indices in the mean wavelength (`wave`) are between the individal channel boundaries
                index = np.where(np.bitwise_and(wave >= wave_low[i], wave <= wave_hi[i]))[0]
                # index = np.where((wave >= wave_low[i])*(wave <= wave_hi[i]))[0] # FINDME: TEST LATER
            else:
                # UVIS: Use all pixels for aperture photometry
                index = range(len(ev.spectra[0,0,0]))
            
            for m in range(ev.n_files):
                '''
                # This is a different way to compute the indices to associate with columns to be summed into 1D spectra
                # Select appropriate orbit-dependent wavelength
                if ev.n_img == (np.max(ev.orbitnum)+1):
                    j = int(ev.orbitnum[m])
                else:
                    j = 0
                #Method 1
                ev.wave.append(np.mean(ev.wavegrid[j][n],axis=0))
                index = np.where(np.bitwise_and(ev.wave[n] >= wave_low, ev.wave[n] <= wave_hi))[0]
                #Method 2
                index = np.where(np.bitwise_and(ev.wave[j][n] >= wave_low, ev.wave[j][n] <= wave_hi))[0]
                '''
                
                # This creates a light curve per NDR
                ev.photflux[n,m]    = np.sum(ev.spectra[m,n,istart:,index],axis=0) # Summing in the 1D spectral plane (m == 1 for WFC3)
                ev.photfluxerr[n,m] = np.sqrt(np.sum(specerr[m,n,istart:,index]**2,axis=0)) # Summing in quadrature the 1D spectral plane (m == 1 for WFC3)
        
        # Save results for individual channels into individual files
        ev.eventname  = ev.eventname2 + '_' + str(int(ev.wave_low*1e3)) + '_' + str(int(ev.wave_hi*1e3))
        # me.saveevent(ev, eventdir + '/d-' + ev.eventname + '-w3', delete=['data_mhdr', 'spectra', 'specerr'])
        
        # saveevent stores everything (that we want) into a pickle
        me.saveevent(ev, eventdir + '/d-' + ev.eventname + '-w3')
        
        # Produce plot
        if isplots == True:
            # 3XYZ: 3: w3 function
            #       X: Plot Number
            #       Y & Z: Spectral Channel Number
            #
            # Normalized Flux vs Time
            plt.figure(3000+i, figsize=(10,8))
            plt.clf() # this clears the frame
            plt.suptitle('Wavelength range: ' + str(wave_low[i]) + '-' + str(wave_hi[i]))
            ax = plt.subplot(111)
            #for n in range(ev.n_spec):
            #plt.subplot(ev.n_spec,1,1)
            #plt.title('Star ' + str(n))
            #igood   = np.where(ev.goodmask[0])[0]
            iscan0  = np.where(ev.scandir == 0)[0]
            iscan1  = np.where(ev.scandir == 1)[0]
            mjd     = np.floor(ev.bjdtdb[0])
            flux0   = np.sum(ev.photflux[0][iscan0],axis=1)/np.sum(ev.photflux[0,[iscan0[-1]]])
            #err  = np.sqrt(1 / np.sum(1/ev.photfluxerr[0]**2,axis=1))/np.sum(ev.photflux[0,-1])
            try:
                err0    = np.sqrt(np.sum(ev.photfluxerr[0][iscan0]**2,axis=1))/np.sum(ev.photflux[0,[iscan0[-1]]])
            except:
                err0    = 0
                # err1    = 0
            plt.errorbar(ev.bjdtdb[iscan0]-mjd, flux0, err0, fmt='bo')
            plt.text(0.05, 0.1, "MAD = "+str(np.round(1e6*np.median(np.abs(np.ediff1d(flux0)))))+" ppm", transform=ax.transAxes, color='b')
            if len(iscan1) > 0:
                flux1   = np.sum(ev.photflux[0][iscan1],axis=1)/np.sum(ev.photflux[0,[iscan0[-1]]])
                err1    = np.sqrt(np.sum(ev.photfluxerr[0][iscan1]**2,axis=1))/np.sum(ev.photflux[0,[iscan1[-1]]])
                plt.errorbar(ev.bjdtdb[iscan1]-mjd, flux1, err1, fmt='ro')
                plt.text(0.05, 0.05, "MAD = "+str(np.round(1e6*np.median(np.abs(np.ediff1d(flux1)))))+" ppm", transform=ax.transAxes, color='r')
            plt.ylabel('Normalized Flux')
            plt.xlabel('Time [MJD + ' + str(mjd) + ']')
            
            plt.subplots_adjust(left=0.10,right=0.95,bottom=0.10,top=0.90,hspace=0.20,wspace=0.3)
            plt.savefig(eventdir + '/figs/' + ev.eventname + '-Fig' + str(3000+i) + '.png')
            #plt.pause(0.1)
            
            if ev.detector == 'IR':
                # Drift: frame number vs drift in the wavelength direction
                plt.figure(3100+i, figsize=(10,8))
                plt.clf()
                for i in range(istart,ev.n_reads-1):
                    plt.subplot(1,np.max((1,ev.n_reads-2)),np.max((1,i)))
                    plt.plot(ev.drift[0,:,i],'.')
                    if i == istart:
                        plt.ylabel('Spectrum Drift')
                    if i == (ev.n_reads-1)/2:
                        plt.xlabel('Frame Number')
                plt.savefig(eventdir + '/figs/' + ev.eventname + '-Fig' + str(3100) + '.png')
        
    if (isplots == True) and (ev.detector == 'IR'):
        # 2D light curve with drift correction
        # Plot frame number vs wavelength with color associated with value
        #   Very cool plot that produces "image" of the entire time series (hopefully corrected)
        plt.figure(3200, figsize=(8,ev.n_files/20.+0.8))
        plt.clf()
        vmin        = 0.97
        vmax        = 1.03
        # istart      = 0
        normspec    = np.mean(ev.spectra[:,0,istart:],axis=1)/np.mean(ev.spectra[-6:,0,istart:],axis=(0,1))
        ediff       = np.zeros(ev.n_files)
        iwmin       = np.where(ev.wave[0][0]>wmin)[0][0]
        iwmax       = np.where(ev.wave[0][0]>wmax)[0][0]
        for i in range(ev.n_files):
            ediff[i]    = 1e6*np.median(np.abs(np.ediff1d(normspec[i,iwmin:iwmax])))
            plt.scatter(ev.wave[0][0], np.zeros(ev.specsize)+i, c=normspec[i], 
                        s=14,linewidths=0,vmin=vmin,vmax=vmax,marker='s',cmap=plt.cm.RdYlBu_r)
        plt.title("MAD = "+str(np.round(np.mean(ediff),0)) + " ppm")
        plt.xlim(wmin,wmax)
        if nchan > 1:
            xticks  = np.round([i for i in np.linspace(wmin, wmax, nchan+1)],3)
            plt.xticks(xticks,xticks)
            plt.vlines(xticks,0,ev.n_files,'k','dashed')
        plt.ylim(0,ev.n_files)
        plt.ylabel('Frame Number')
        plt.xlabel('Wavelength ($\mu m$)')
        plt.colorbar()
        plt.tight_layout()
        plt.savefig(eventdir+'/figs/fig3200-2D_LC.png')

#lcGMOS(isplots=True)

'''
flux1 = np.sum(ev.photflux[0:6],axis=0)
flux2 = np.sum(ev.photflux[6:12],axis=0)
flux3 = np.sum(ev.photflux[12:18],axis=0)
diffflux = flux2/(flux1+flux3)
normflux = diffflux/np.median(diffflux)
plt.figure(1)
plt.clf()
plt.suptitle('WASP-12 Light Curve')
plt.plot(ev.bjd_tdb-ev.bjd_tdb[0], normflux, 'o')
plt.xlabel('Time (Days)')
plt.ylabel('Normalized Flux')
plt.savefig('figs/' + ev.eventname + '-LC.png')
'''
