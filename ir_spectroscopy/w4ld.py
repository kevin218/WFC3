#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy as sp
import scipy.interpolate as spi
from astropy.io import fits
from ..lib import manageevent as me
from ..lib import getldparam as gld
from ..lib import limbdark

from ..lib import integrate

# Calculate limb-darkening coefficients using Phoenix data
def limbDarkening(filename, wavelow, wavehi, specwave=None, spectrum=None, n_param=4, n_plot=False, stellarmodel='phoenix'):
    '''

    '''
    # Kurucz stellar models
    if stellarmodel == 'kurucz':
        print("Using Kurucz stellar models.")
        spit = np.array((specwave, spectrum))
        coeffs, mu, sum_int = gld.getlimbparam(filename, 0, n_param, lamfact=0.1, spit=spit, returnmu=True)
        #print(coeffs)
        if n_param == 1:
            model = limbdark.linear(mu, coeffs)
        if n_param == 2:
            model = limbdark.quad(mu, coeffs)
        if n_param == 4:
            model = limbdark.nl(mu, coeffs)
        if n_plot != False:
            plt.figure(n_plot)
            plt.clf()
            plt.plot(mu, sum_int, 'go', label='Intensities')
            plt.plot(mu, model, 'r-', label='Best-Fit Model')
            plt.xlabel('mu')
            plt.ylabel('Normalized Intensity')
            plt.legend(loc='lower right')

    # Phoenix stellar models
    elif stellarmodel == 'phoenix':
        print("Using Phoenix stellar models.")
        # Read file
        hdulist = fits.open(filename)
        lddata  = hdulist[1].data

        # Read data
        wave      = lddata['wavelength'] / 1e4
        igood     = np.where(np.bitwise_and(wave >= wavelow, wave <= wavehi))
        wave      = wave[igood]
        intensity = lddata['intensity'][igood].T * wave
        flux      = lddata['flux'][igood]
        mu        = lddata['intensity'][0]
        n_mu, n_weights = intensity.shape

        if wave != None and spectrum != None:
            # Weight intensities according to the spectrum
            spline      = spi.UnivariateSpline(specwave, spectrum)
            newspectrum = spline(wave)
            intbottom   = integrate.integrate(specwave, spectrum)
            sum_int     = np.zeros(n_mu)
            for i in range(n_mu):
                inttop     = integrate.integrate(wave, intensity[i]/intensity[-1]*newspectrum)
                sum_int[i] = inttop / intbottom
        else:
            # Sum intensities at each mu
            sum_int  = np.sum(intensity, axis=1)

        sum_int /= sum_int.max()

        if n_plot != False:
            plt.figure(n_plot)
            plt.clf()
            plt.plot(mu, sum_int, 'bo', label='Original Intensities')

        # Set stellar limb at 1/e of max intensity
        spline   = spi.UnivariateSpline(sum_int, mu, k=1, s=0)
        mu_zero  = spline(1./np.e)
        scale    = 1./np.cos(np.pi/2 - np.arccos(mu_zero))
        isstar   = np.where(mu >= mu_zero)[0]

        # Recalculate angles and intensities
        mu       = np.sqrt(1. - (1. - mu[isstar]**2) * scale**2)
        n_mu     = mu.size
        sum_int  = sum_int[isstar]

        #Trim mu < 0.1
        imu     = np.where(mu >= 0.1)[0]
        mu      = mu[imu]
        sum_int = sum_int[imu]

        if n_plot != False:
            plt.figure(n_plot)
            plt.plot(mu, sum_int, 'go', label='Rescaled Intensities')

        # Non-Linear Limb Darkening
        # I(mu) / I(1) = 1 - a1(1 - mu^[1/2]) - a2(1 - mu)
        #                  - a3(1 - mu^[3/2]) - a4(1 - mu^2)
        # Therefore:
        # 1 - I(mu) / I(1) =   a1(1 - mu^[1/2]) + a2(1 - mu)
        #                    + a3(1 - mu^[3/2]) + a4(1 - mu^2)
        #
        # Quadratic Limb Darkening
        # I(mu) / I(1) = 1 - a1(1 - mu) - a2(1 - mu)^2

        # Use least sqaures to find parameters (a1, etc.) for best
        # model fit using the above equation
        A       = np.zeros((n_mu, n_param))

        #Linear
        if   n_param == 1:
            A[:, 0] = (1 - mu)
            coeffs = np.linalg.lstsq(A, 1 - sum_int)[0]
            model = limbdark.linear(mu, coeffs)
        #Quadratic
        if   n_param == 2:
            A[:, 0] = (1 - mu)
            A[:, 1] = (1 - mu)**2
            coeffs = np.linalg.lstsq(A, 1 - sum_int)[0]
            model = limbdark.quad(mu, coeffs)
        #Non Linear
        elif n_param == 4:
            A[:, 0] = (1 - mu**0.5)
            A[:, 1] = (1 - mu     )
            A[:, 2] = (1 - mu**1.5)
            A[:, 3] = (1 - mu**2. )
            coeffs = np.linalg.lstsq(A, 1 - sum_int)[0]
            model = limbdark.nl(mu, coeffs)
        else:
            print('WARNING: ' + str(n_param) + ' parameter models are not supported.')
            coeffs = np.zeros(n_param)
            model  = np.zeros(n_mu)

        if n_plot != False:
            plt.figure(n_plot)
            plt.plot(mu, model, 'r-', label='Best-Fit Model')
            plt.xlabel('mu')
            plt.ylabel('Normalized Intensity')
            plt.legend(loc='lower right')

    else:
        print("Unrecognized stellar models.")
        coeffs = None

    return coeffs


def ld_driver(eventname, eventdir, wave_low=None, wave_hi=None, n_param=4, isplots=False, stellarmodel='phoenix'):
    '''

    '''
    # Load saved data
    print("Loading saved data...")
    ev              = me.loadevent(eventdir + '/d-' + eventname + '-w3')
    aux             = me.loadevent(eventdir + '/d-' + ev.eventname2 + '-data')
    ev.spectra      = aux.spectra

    '''
    #FINDME: HACK
    ev.file_med = ev.loc_ld + 'kelt11_med.txt'
    ev.file_med = ev.loc_ld + 'lte6250-4.38+0.2a+0.0CMg-0.1.BT-dusty-giant-2011.cifist.He.irf.fits'
    ev.file_low = ev.loc_ld + 'lte6100-4.38+0.2a+0.0CMg-0.1.BT-dusty-giant-2011.cifist.He.irf.fits'
    ev.file_hi  = ev.loc_ld + 'lte6400-4.38+0.2a+0.0CMg-0.1.BT-dusty-giant-2011.cifist.He.irf.fits'
    #print(ev.file_med)
    '''
    n   = 0
    m   = 0
    ilo = np.where(ev.wave[n][m] > wave_low)[0][0]
    ihi = np.where(ev.wave[n][m] < wave_hi)[0][-1]+1

    #
    print("Computing limb-darkening coefficients...")
    specwave = ev.wave[n][m][ilo:ihi]*1e4   #Angstroms
    #iwave    = np.argsort(specwave)
    #specwave = specwave[iwave]*1e4  #Angstroms
    wavelow  = specwave[0]
    wavehi   = specwave[-1]
    spectrum = np.sum(ev.spectra[n,:,ilo:ihi],axis=0)
    if isplots:
        # Optimal
        ev.ldcoeffs     = limbDarkening(ev.file_med, wavelow, wavehi, specwave=specwave, spectrum=spectrum, n_param=n_param, n_plot=4000, stellarmodel=stellarmodel)
        plt.title(str(n_param) + ' parameter model, optimal fit')
        plt.savefig(eventdir + '/figs/' + ev.eventname + '-Fig' + str(4000) + str(n_param) + '.png')
        try:
            # Low
            ev.ldcoeffs_low = limbDarkening(ev.file_low, wavelow, wavehi, specwave=specwave, spectrum=spectrum, n_param=n_param, n_plot=4001, stellarmodel=stellarmodel)
            plt.title(str(n_param) + ' parameter model, low fit')
            plt.savefig(eventdir + '/figs/' + ev.eventname + '-Fig' + str(4001) + str(n_param) + '.png')
            # Hi
            ev.ldcoeffs_hi  = limbDarkening(ev.file_hi, wavelow, wavehi, specwave=specwave, spectrum=spectrum, n_param=n_param, n_plot=4002, stellarmodel=stellarmodel)
            plt.title(str(n_param) + ' parameter model, hi fit')
            plt.savefig(eventdir + '/figs/' + ev.eventname + '-Fig' + str(4002) + str(n_param) + '.png')
        except:
            pass
    else:
        ev.ldcoeffs     = limbDarkening(ev.file_med, wavelow, wavehi, specwave=specwave, spectrum=spectrum, n_param=n_param, n_plot=False, stellarmodel=stellarmodel)
        try:
            ev.ldcoeffs_low = limbDarkening(ev.file_low, wavelow, wavehi, specwave=specwave, spectrum=spectrum, n_param=n_param, n_plot=False, stellarmodel=stellarmodel)
            ev.ldcoeffs_hi  = limbDarkening(ev.file_hi, wavelow, wavehi, specwave=specwave, spectrum=spectrum, n_param=n_param, n_plot=False, stellarmodel=stellarmodel)
        except:
            pass
    print(eventname, ev.ldcoeffs)

    # Save results
    print('Saving results...')
    me.saveevent(ev, eventdir + '/d-' + ev.eventname + '-w4', delete=['spectra'])

    return

#ld_driver(True)
