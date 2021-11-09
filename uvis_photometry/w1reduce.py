#! /usr/bin/env python

import numpy as np
import scipy.interpolate as spi
from astropy.io import fits
import matplotlib.pyplot as plt
import multiprocessing as mp
import time, sys, os, shutil
from ..lib import manageevent as me
from ..lib import sort_nicely as sn
from ..lib import hst_scan as hst
from ..lib import optspex, centroid, suntimecorr, utc_tt

#NOTE: hst_scan works for stare mode when using ima files
#import hst
import scipy.ndimage.interpolation as spni

try:
  basestring
except NameError:
  basestring = str
 
def reduceWFC3(eventname, isplots=False):
    '''
    Reduces data images and calculated optimal spectra.
    
    Parameters
    ----------
    isplots     : Set True to produce plots
    
    Returns
    -------
    None
    
    Remarks
    -------
    Requires eventname_params file to intialize event object
    
    History
    -------
    Written by Kevin Stevenson      June 2012
    
    '''
    evpname = eventname + '_params'
    exec("import {0} as evp".format(evpname), globals())
    # exec('import ' + evpname + ' as evp')
    # reload(evp)
    
    t0 = time.time()
    # Initialize event object
    # All parameters are specified in this file
    ev = evp.event_init()
    try:
        aux = evp.aux_init()
    except:
        print("Need to update event file to include auxiliary object.")
        return
        
    # Create directories
    if not os.path.exists(ev.eventdir):
        os.makedirs(ev.eventdir)
    if not os.path.exists(ev.eventdir+"/figs"):
        os.makedirs(ev.eventdir+"/figs")
    
    # Copy ev_params file
    shutil.copyfile(evpname + '.py', ev.eventdir+'/'+evpname+'.py')
    
    # Load frame filenames
    # Flat
    '''
    ev.flat_list = []
    if ev.flatfile == None:
        for i in range(ev.flatstart,ev.flatend):
            ev.flat_list.append(ev.loc_cal + ev.filebase + str(i).zfill(4) + ".fits")
    else:
        handle = open(ev.flatfile)
        for line in handle:
            ev.flat_list.append(ev.loc_cal + line)
        handle.close()
    # Bias
    ev.bias_list = []
    if ev.biasfile == None:
        for i in range(ev.biasstart,ev.biasend):
            ev.bias_list.append(ev.loc_cal + ev.filebase + str(i).zfill(4) + ".fits")
    else:
        handle = open(ev.biasfile)
        for line in handle:
            ev.bias_list.append(ev.loc_cal + line)
        handle.close()
    # wavelenth
    ev.wave_list = []
    if ev.wavefile == None:
        for i in range(ev.wavestart,ev.waveend):
            ev.wave_list.append(ev.loc_cal + ev.filebase + str(i).zfill(4) + ".fits")
    else:
        handle = open(ev.wavefile)
        for line in handle:
            ev.wave_list.append(ev.loc_cal + line)
        handle.close()
    '''
    # Object
    ev.obj_list = []   #Do not rename ev.obj_list!
    if ev.objfile == None:
        #Retrieve files within specified range
        for i in range(ev.objstart,ev.objend):
            ev.obj_list.append(ev.loc_sci + ev.filebase + str(i).zfill(4) + ".fits")
    elif ev.objfile == 'all':
        #Retrieve all files from science directory
        for fname in os.listdir(ev.loc_sci):
            ev.obj_list.append(ev.loc_sci +'/'+ fname)
        ev.obj_list = sn.sort_nicely(ev.obj_list)
    else:
        #Retrieve filenames from list
        handle = open(ev.objfile)
        for line in handle:
            ev.obj_list.append(ev.loc_sci + line)
        handle.close()
    ev.n_files = len(ev.obj_list)
    
    #Determine image size and filter/grism
    hdulist         = fits.open(ev.obj_list[0])
    nx              = hdulist['SCI',1].header['NAXIS1']
    ny              = hdulist['SCI',1].header['NAXIS2']
    ev.grism        = hdulist[0].header['FILTER']
    ev.detector     = hdulist[0].header['DETECTOR']
    ev.flatoffset   = [[-1*hdulist['SCI',1].header['LTV2'], -1*hdulist['SCI',1].header['LTV1']]]
    hdulist.close()
    
    # Calculate centroid of direct image(s)
    ev.img_list = []
    if isinstance(ev.directfile, basestring) and ev.directfile.endswith('.fits'):
        ev.img_list.append(ev.loc_cal + ev.directfile)
    else:
        #Retrieve filenames from list
        handle = open(ev.directfile)
        for line in handle:
            ev.img_list.append(ev.loc_cal + line)
        handle.close()
    ev.n_img = len(ev.img_list)
    ev.centroid, ev.directim = hst.imageCentroid(ev.img_list, ev.centroidguess, ev.centroidtrim, ny, ev.obj_list[0])
    
    # Calculate theoretical centroids along spatial scan direction
    ev.centroids = []
    for j in range(ev.n_img):
        ev.centroids.append([])
        for i in range(ev.n_spec):
            #FINDME: Can't assume that scan direction is only in y direction (no x component)
            ev.centroids[j].append([np.zeros(ev.ywindow[i][1]-ev.ywindow[i][0])+ev.centroid[j][0],ev.centroid[j][1]])
            #ev.centroids[j].append([np.arange(ev.ywindow[i][0]-5,ev.ywindow[i][1]),ev.centroid[j][1]])
            #Apply small shift in y direction to allow griddata to work
            #griddata requires that three points can't be perfectly aligned
            #ev.centroids[j][i][0] = ev.centroids[j][i][0] + np.random.normal(0,1e-7,len(ev.centroids[j][i][0]))
    
    # Calculate trace
    #refpix = [152,60]
    #yoff, xoff = refpix - ev.centroid
    print("Calculating 2D trace and wavelength assuming " + ev.grism + " filter/grism...")
    ev.xrange   = []
    ev.yrange   = []
    ev.yrange5  = []
    for i in range(ev.n_spec):
        ev.xrange.append(np.arange(ev.xwindow[i][0],ev.xwindow[i][1]))
        ev.yrange.append(np.arange(ev.ywindow[i][0],ev.ywindow[i][1]))
        ev.yrange5.append(np.arange(ev.ywindow[i][0]-5,ev.ywindow[i][1]))
    ev.trace2d  = []
    ev.wave2d   = []
    for j in range(ev.n_img):
        ev.trace2d.append([])
        ev.wave2d.append([])
        for i in range(ev.n_spec):
            ev.trace2d[j].append(hst.calcTrace(ev.xrange[i], ev.centroids[j][i], ev.grism))
            # print(ev.xrange[i], ev.centroids[j][i], ev.grism)
            ev.wave2d[j].append(hst.calibrateLambda(ev.xrange[i], ev.centroids[j][i], ev.grism)/1e4)     #wavelength in microns
    
    if ev.detector == 'IR':
        print("Calculating slit shift values...")
        i = 0   # Use first spectrum
        j = -1  # Use last image
        spectrum    = fits.getdata(ev.obj_list[j])
        ev.slitshift, ev.shift_values, ev.yfit = hst.calc_slitshift2(spectrum, ev.xrange[i], ev.ywindow[i], ev.xwindow[i])
        ev.wavegrid  = ev.wave2d
        ev.wave = []
        for j in range(ev.n_img):
            ev.wave.append([])
            for i in range(ev.n_spec):
                # print(type(ev.wavegrid))
                # print(type(ev.wave))
                # print(np.shape(ev.wavegrid))
                # print(ev.wave[j])
                # print(ev.wavegrid[j][i])
                ev.wave[j].append(np.mean(ev.wavegrid[j][i],axis=0))
    else:
        # FINDME: Assume no slitshift
        ev.yfit         = range(ev.ywindow[0][1] - ev.ywindow[0][0])
        ev.slitshift    = np.zeros(ev.ywindow[0][1] - ev.ywindow[0][0])
        ev.shift_values = np.zeros(len(ev.yfit))
    #********************
    """
    hdu = fits.open(ev.obj_list[0])
    frames = hdu['sci',14].data
    plt.figure(1)
    plt.clf()
    plt.imshow(frames,aspect='auto',origin='lower',vmin=0)
    plt.plot(ev.xrange[0],ev.trace[0],'k-',lw=3)
    plt.plot([ev.centroid[1]],[ev.centroid[0]],'r+')
    plt.ylim(148,250)
    plt.xlim(20,200)
    
    plt.figure(5)
    plt.clf()
    #plt.imshow(ev.wavegrid[0],origin='lower',aspect='auto')
    for i in [10,49,99,149]:
        plt.plot(ev.xrange[0],ev.wavegrid[0][i],'-')
    
    plt.figure(2)
    plt.clf()
    plt.plot(frames[:,100], 'r-o')
    plt.plot(frames[:,150], 'b-o')
    plt.vlines(ev.trace[0][100-50], 0, 5000, 'r')
    plt.vlines(ev.trace[0][150-50], 0, 5000, 'b')
    plt.xlim(148,250)
    
    plt.figure(3)
    plt.clf()
    plt.imshow(ev.directim,aspect='auto',origin='lower',vmin=0, vmax=1000)
    plt.plot([ev.centroid[1]],[ev.centroid[0]],'r+')
    plt.xlim(10,40)
    plt.ylim(135,165)
    
    plt.figure(4)
    plt.clf()
    plt.plot(ev.wave[0], frames[0,162,50:200], '-')
    #plt.plot(wave, frames[0,154,50:200], '-')
    #plt.plot(wave, frames[0,155,50:200], '-')
    #plt.plot(wave, frames[0,156,50:200], '-')
    #plt.plot(wave, frames[0,157,50:200], '-')
    """
    #********************
    """
    # Make list of master bias frames
    # bias_master must be a list to account for frames of varying size
    if len(ev.bias_list) > 0:
        print('Loading bias frames...')
        bias_frames, biaserr_frames, bias_hdr, bias_mhdr = hst.read(ev.bias_list)
        ev.bias_master = []
        for n in range(ev.n_spec):
            ev.bias_master.append(np.median(bias_frames[n], axis=0))
    else:
        print('No bias frames found.')
        #bias_frames  = np.zeros((ev.n_spec,1,ny,nx))
        ev.bias_master = np.zeros((ny,nx))
        bias_hdr     = None
        bias_mhdr    = None
    """
    print('No bias frames found.')
    ev.bias_master = np.zeros((ny,nx))
    bias_hdr     = None
    bias_mhdr    = None
    
    # Make list of master flat field frames
    if ev.flatfile == None:
        print('No flat frames found.')
        flat_master = np.ones((ev.n_img,ev.n_spec,ny,nx))
        mask        = np.ones((ev.n_img,ev.n_spec,ny,nx))
        flat_hdr    = None
        flat_mhdr   = None
    else:
        print('Loading flat frames...')
        #ev.flatfile = ev.loc_cal+ev.flatfile
        #wave = np.mean(ev.wave,axis=1)
        flat_master = []
        mask = []
        for j in range(ev.n_img):
            flat, tempmask = hst.makeflats(ev.flatfile, ev.wavegrid[j], ev.xwindow, ev.ywindow, ev.flatoffset, ev.n_spec, ny, nx, sigma=ev.flatsigma)
            flat_master.append(flat)
            mask.append(tempmask)
    
    # Manually mask regions [specnum, colstart, colend]
    if hasattr(ev, 'manmask'):
        print("\rMasking manually identified bad pixels.")
        for j in range(ev.n_img):
            for i in range(len(ev.manmask)):
                ind, colstart, colend, rowstart, rowend = ev.manmask[i]
                n = ind % ev.n_spec
                mask[j][n][rowstart:rowend,colstart:colend] = 0 #ev.window[:,ind][0]:ev.window[:,ind][1]
    
    if isplots:
        # Plot normalized flat fields
        plt.figure(1000, figsize=(12,8))
        plt.clf()
        plt.suptitle('Master Flat Frames')
        for i in range(ev.n_spec):
            for j in range(ev.n_img):
                #plt.subplot(ev.n_spec,ev.n_img,i*ev.n_img+j+1)
                plt.subplot(2,np.ceil(ev.n_img/2.),i*ev.n_img+j+1)
                plt.title(str(j) +','+ str(i))
                plt.imshow(flat_master[j][i], origin='lower')
        plt.savefig(ev.eventdir + '/figs/fig1000-Flats.png')
        # Plot masks
        plt.figure(1001, figsize=(12,8))
        plt.clf()
        plt.suptitle('Mask Frames')
        for i in range(ev.n_spec):
            for j in range(ev.n_img):
                #plt.subplot(ev.n_spec,ev.n_img,i*ev.n_img+j+1)
                plt.subplot(2,np.ceil(ev.n_img/2.),i*ev.n_img+j+1)
                plt.title(str(j) +','+ str(i))
                plt.imshow(mask[j][i], origin='lower')
        plt.savefig(ev.eventdir + '/figs/fig1001-Masks.png')
        """
        # Plot master bias frames
        plt.figure(1001, figsize=(12,8))
        plt.clf()
        plt.suptitle('Master Bias Frames')
        for n in range(ev.n_spec):
            plt.subplot(1,ev.n_spec,n+1)
            #plt.title(str(n+1))
            plt.imshow(ev.bias_master, origin='lower')
        plt.savefig(ev.eventdir + '/figs/fig1001-Bias.png')
        # Plot Slit shift
        plt.figure(1004, figsize=(12,8))
        plt.clf()
        plt.suptitle('Model Slit Tilts/Shifts')
        for i in range(ev.n_spec):
            for j in range(ev.n_img):
                plt.subplot(ev.n_spec,ev.n_img,i*ev.n_img+j+1)
                plt.plot(ev.shift_values[j][i], range(np.size(ev.shift_values[j][i])), '.')
                plt.plot(ev.slitshift[j][i], range(np.size(ev.slitshift[j][i])), 'r-', lw=2)
        plt.savefig(ev.eventdir + '/figs/fig1004-SlitTilt.png')
        plt.pause(0.1)
        """
        
        # Plot Slit shift
        plt.figure(1004, figsize=(12,8))
        plt.clf()
        plt.suptitle('Model Slit Tilts/Shifts')
        plt.plot(ev.shift_values, ev.yfit, '.')
        plt.plot(ev.slitshift, range(ev.ywindow[0][0],ev.ywindow[0][1]), 'r-', lw=2)
        plt.xlim(-1,1)
        plt.savefig(ev.eventdir + '/figs/fig1004-SlitTilt.png')
        plt.pause(0.1)
    
    # Write spectrum to file
    global spectra, specerr, specbg, fracMaskReg, fracMaskOpt, scanHeight
    def writeSpectrum(arg):
        '''
        
        '''
        spectrum, specstd, stdbg, fmreg, fmopt, scanH, m, n = arg
        # Assign to array of spectra and uncertainties
        spectra[m,n]    = spectrum
        specerr[m,n]    = specstd
        specbg[m,n]     = stdbg
        fracMaskReg[m,n] = fmreg
        fracMaskOpt[m,n] = fmopt
        scanHeight[m,n] = scanH
        return

    # Determine if we are using IMA or FLT files
    # FLT files subtract first from last, 2 reads
    if ev.obj_list[0].endswith('flt.fits'):
        ev.n_reads  = 2
    else:
        hdulist     = fits.open(ev.obj_list[0])
        ev.n_reads  = hdulist['SCI',1].header['SAMPNUM']
        hdulist.close()
    # Create arrays
    ev.specsize = ev.xwindow[0][1]-ev.xwindow[0][0]
    spectra     = np.zeros((ev.n_files, ev.n_spec, ev.n_reads-1, ev.specsize))
    specerr     = np.zeros((ev.n_files, ev.n_spec, ev.n_reads-1, ev.specsize))
    specbg      = np.zeros((ev.n_files, ev.n_spec, ev.n_reads-1, ev.specsize))
    fracMaskReg = np.zeros((ev.n_files, ev.n_spec, ev.n_reads-1))
    fracMaskOpt = np.zeros((ev.n_files, ev.n_spec, ev.n_reads-1))
    scanHeight  = np.zeros((ev.n_files, ev.n_spec))
    #frames, err, data_hdr, data_mhdr = hst.read(obj_list)
    
    # Record JD and exposure times
    print('Reading headers, recording JD and exposure times...')
    ev.jd       = np.zeros(ev.n_files)
    ev.exptime  = np.zeros(ev.n_files)
    data_mhdr   = []
    data_hdr   = []
    
    for m in range(ev.n_files):
        data, err, hdr, mhdr = hst.read(ev.obj_list[m])
        data_mhdr.append(mhdr[0])
        data_hdr.append(hdr[0])
        ev.jd[m]       = 2400000.5 + 0.5*(data_mhdr[m]['EXPSTART'] + data_mhdr[m]['EXPEND'])
        ev.exptime[m]  = data_mhdr[m]['EXPTIME']
    del(data,err)
    
    ev.readNoise = np.mean((data_mhdr[0]['READNSEA'],
                            data_mhdr[0]['READNSEB'],
                            data_mhdr[0]['READNSEC'],
                            data_mhdr[0]['READNSED']))
    print('Read noise: ' + str(ev.readNoise))
    print('Gain: ' + str(ev.gain))
    #ev.v0 = (ev.readNoise/ev.gain)**2     #Units of ADU
    ev.v0 = ev.readNoise**2                #Units of electrons
    
    # Group frames into frame, batch, and orbit number
    ev.framenum, ev.batchnum, ev.orbitnum = hst.groupFrames(ev.jd)
    
    ev.ra       = data_mhdr[0]['RA_TARG']*np.pi/180
    ev.dec      = data_mhdr[0]['DEC_TARG']*np.pi/180
    if ev.horizonsfile != None:
        # Apply light-time correction, convert to BJD_TDB
        # Horizons file created for HST around time of observations
        print("Converting times to BJD_TDB...")
        ev.bjd_corr = suntimecorr.suntimecorr(ev.ra, ev.dec, ev.jd, ev.horizonsfile)
        bjdutc      = ev.jd + ev.bjd_corr/86400.
        ev.bjdtdb   = utc_tt.utc_tt(bjdutc)
        print('BJD_corr range: ' + str(ev.bjd_corr[0]) + ', ' + str(ev.bjd_corr[-1]))
    else:
        print("No Horizons file found.")
        ev.bjdtdb   = ev.jd
    
    #mask = np.ones((ev.n_spec,ev.n_files,ny,nx))
    print('Processing data frames...')
    print('Total # of files: '+str(ev.n_files))
    if ev.ncpu == 1:
        # Only 1 CPU
        for m in range(ev.n_files):
            #Select appropriate flat, mask, and slitshift
            if ev.n_img == (np.max(ev.orbitnum)+1):
                j = int(ev.orbitnum[m])
            else:
                j = 0
            for n in range(ev.n_spec):
                writeSpectrum(hst.calcSpectrum(ev.obj_list[m], mask[j][n], ev.bias_master, flat_master[j][n], ev.slitshift, ev.xwindow[n], ev.ywindow[n], ev.gain, ev.v0, ev.spec_width, ev.fitbghw, m, n, diffthresh=ev.diffthresh, p3thresh=ev.p3thresh, p5thresh=ev.p5thresh, p7thresh=ev.p7thresh, fittype=ev.fittype, window_len=ev.window_len, deg=ev.deg, expand=ev.expand, isplots=isplots, eventdir=ev.eventdir, bgdeg=ev.bgdeg))
    else:
        # Multiple CPUs
        pool = mp.Pool(ev.ncpu)
        for m in range(ev.n_files):
            #Select appropriate flat, mask, and slitshift
            if ev.n_img == (np.max(ev.orbitnum)+1):
                j = int(ev.orbitnum[m])
            else:
                j = 0
            for n in range(ev.n_spec):
                res = pool.apply_async(hst.calcSpectrum, args=(ev.obj_list[m], mask[j][n], ev.bias_master, flat_master[j][n], ev.slitshift, ev.xwindow[n], ev.ywindow[n], ev.gain, ev.v0, ev.spec_width, ev.fitbghw, m, n, ev.diffthresh, ev.p3thresh, ev.p5thresh, ev.p7thresh, ev.fittype, ev.window_len, ev.deg, ev.expand, False, ev.eventdir, ev.bgdeg), callback=writeSpectrum)
        
        pool.close()
        pool.join()
        res.wait()

    # Calculate total time
    total = (time.time() - t0)/60.
    print('\nTotal time (min): ' + str(np.round(total,2)))
    
    # Save results
    print('Saving results...')
    ev.scanHeight   = scanHeight
    aux.spectra     = spectra
    aux.specerr     = specerr
    aux.specbg      = specbg
    aux.fracMaskReg = fracMaskReg
    aux.fracMaskOpt = fracMaskOpt
    aux.data_hdr    = data_hdr
    aux.data_mhdr   = data_mhdr
    aux.mask        = mask
    aux.bias_mhdr   = bias_mhdr
    aux.flat_master = flat_master
    #np.savez_compressed(ev.eventdir + '/d-'+ev.eventname+'-data', spectra=spectra, specerr=specerr, specbg=specbg, data_hdr=data_hdr, data_mhdr=data_mhdr, mask=mask, bias_mhdr=bias_mhdr, flat_master=flat_master) #frames=frames, 
    me.saveevent(aux, ev.eventdir + '/d-' + ev.eventname + '-data')
    me.saveevent( ev, ev.eventdir + '/d-' + ev.eventname + '-w1')    #, delete=['flat_mhdr'])

    if ev.detector == 'IR':
        # 2D light curve without drift correction
        plt.figure(1012, figsize=(8,ev.n_files/20.+0.8))
        plt.clf()
        vmin        = 0.96
        vmax        = 1.04
        istart      = 1
        normspec    = np.mean(spectra[:,0,istart:],axis=1)/np.mean(spectra[ev.inormspec[0]:ev.inormspec[1],0,istart:],axis=(0,1))
        for i in range(ev.n_files):
            plt.scatter(ev.wave[0][0], np.zeros(ev.specsize)+i, c=normspec[i], 
                        s=14,linewidths=0,vmin=vmin,vmax=vmax,marker='s',cmap=plt.cm.RdYlBu_r)
        plt.xlim(1.1,1.7)
        plt.ylim(0,ev.n_files)
        plt.ylabel('Frame Number')
        plt.xlabel('Wavelength ($\mu m$)')
        plt.colorbar()
        plt.tight_layout()
        plt.savefig(ev.eventdir+'/figs/fig1012-2D_LC.png')
    else:
        # 2D light curve without drift correction
        plt.figure(1012, figsize=(8,ev.n_files/20.+0.8))
        plt.clf()
        vmin        = 0.96
        vmax        = 1.04
        normspec    = spectra[:,0,0]/np.mean(spectra[ev.inormspec[0]:ev.inormspec[1],0,0],axis=0)
        for i in range(ev.n_files):
            plt.scatter(range(ev.specsize), np.zeros(ev.specsize)+i, c=normspec[i], 
                        s=14,linewidths=0,vmin=vmin,vmax=vmax,marker='s',cmap=plt.cm.RdYlBu_r)
        plt.xlim(0,ev.specsize)
        plt.ylim(0,ev.n_files)
        plt.ylabel('Frame Number')
        plt.xlabel('Wavelength ($\mu m$)')
        plt.colorbar()
        plt.tight_layout()
        plt.savefig(ev.eventdir+'/figs/fig1012-2D_LC.png')
    
    return

'''
handle  = np.load('wa012-reduced.npz')
spectra = handle['spectra']
specerr = handle['specerr']
handle.close()

n_files, n_spec, specsize = spectra.shape
for m in range(n_files):
    print(m,n_files)
    plt.figure(1000+m, figsize=(20,8))
    plt.clf()
    for n in range(n_spec):
        plt.subplot(3,6,n+1)
        plt.errorbar(range(3,specsize-3), spectra[m,n,3:-3], specerr[m,n,3:-3], fmt='-')
        plt.xlim(0,520)
        plt.subplots_adjust(left=0.05,right=0.98,bottom=0.10,top=0.95,hspace=0.20,wspace=0.3)
    plt.savefig('figs/' + obj + '-Fig' + str(1000+m) + '.png')
    plt.close(1000+m)

'''

