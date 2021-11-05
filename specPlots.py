
import os, sys, copy
import manageevent as me
from run import *

#filestr = '/Users/kristin/Documents/STScI/KBS_WFC3-IR_Algorithm/2017-07-15-w1_spec_width_20/fitbghw_20/divideWLC_09272017/2017-09-28_16-34-test/'
#filestr = '/Users/kristin/Documents/STScI/KBS_WFC3-IR_Algorithm/2017-07-15-w1_spec_width_20/fitbghw_20/divideWLC_09072017/2017-09-18_08-26-test/'
filestr = '/Users/kristin/Documents/STScI/KBS_WFC3-IR_Algorithm/2017-07-15-w1_spec_width_20/fitbghw_20/spectroscopicLC_09272017/2017-09-28_16-28-test/'
ev = w6Restore(filedir = filestr)
'''
ev1 = w6Restore(filedir = filestr)
ev2 = w6Restore(filedir = filestr)
ev3 = w6Restore(filedir = filestr)
ev4 = w6Restore(filedir = filestr)
ev5 = w6Restore(filedir = filestr)
ev6 = w6Restore(filedir = filestr)
ev7 = w6Restore(filedir = filestr)
ev8 = w6Restore(filedir = filestr)
ev9 = w6Restore(filedir = filestr)

ev          = [ev0[0],ev1[0],ev2[0],ev3[0],ev4[0],ev5[0],ev6[0],ev7[0],ev8[0],ev9[0]]
'''
evall       = [ev[0],ev[1],ev[2],ev[3],ev[4],ev[5],ev[6],ev[7],ev[8],ev[9]]
#evall       = [ev[0][0],ev[1][0],ev[2][0],ev[3][0],ev[4][0],ev[5][0],ev[6][0],ev[7][0],ev[8][0],ev[9][0]]
#evall       = [ev0,ev1,ev2,ev3,ev4,ev5,ev6,ev7,ev8,ev9]

temp = evall[0]
print('Printing ev')
print(ev[0])
#nobs        = len(ev0)
numchan     = len(ev)
wave        = np.zeros(numchan)
waveerr     = np.zeros(numchan)
depth       = np.zeros(numchan)
deptherr    = np.zeros(numchan)
wavelow     = np.zeros(numchan)
wavehi      = np.zeros(numchan)
linearLD    = np.zeros(numchan)
for i in range(numchan):
    wave[i]     = (ev[i].wave_hi+ev[i].wave_low)/2.
    waveerr[i]  = (ev[i].wave_hi-ev[i].wave_low)/2.
    depth[i]    = 100*ev[i].fit[0].bestp[1]**2
    deptherr[i] = 200*ev[i].fit[0].bestp[1]*ev[i].fit[0].medianp[1,1]
    wavelow[i]  = ev[i].wave_low
    wavehi[i]   = ev[i].wave_hi
    linearLD[i] = ev[i].fit[0].bestp[6]
    
#Plot spectrum
plt.figure(1, figsize=(10,6))
plt.clf()
plt.title('WASP-79b', size=16)
plt.errorbar(wave, depth, deptherr, waveerr, fmt='bo', ecolor='k')
plt.ylabel('Transit Depth (%)', size=14)
plt.xlabel('Wavelength ($\mu m$)', size=14)
#plt.ylim(0,200)
plt.subplots_adjust(0.09,0.09,0.98,0.95)

#plt.savefig('WASP79b-WFC3-TrSpec.png')
#plt.savefig('WASP79b-WFC3-TrSpec.pdf')

#Plot spectrum
plt.figure(2, figsize=(10,6))
plt.clf()
plt.title('WASP-79b', size=16)
plt.errorbar(wave, linearLD, deptherr, waveerr, fmt='bo', ecolor='k')
plt.ylabel('Linear LD Parameter', size=14)
plt.xlabel('Wavelength ($\mu m$)', size=14)
#plt.ylim(0,200)
plt.subplots_adjust(0.09,0.09,0.98,0.95)

"""
#Compute reduced Chi-squared for binned data (for verification)
redchi2 = np.zeros(numchan)
for i in range(numchan):
    for j in range(nobs):
        redchi2[i] += evall[i][j].fit[0].redchisq/nobs
    print(i,np.round(redchi2[i],2))



# Combine residuals for RMS vs Binsize plot
res     = []
phase   = []
for i in range(numchan):
    res.append(np.array([]))
    phase.append(np.array([]))
    for j in range(nobs):
        res[i]  = np.concatenate((    res[i], evall[i][j].fit[0].normresiduals))
        phase[i]= np.concatenate((  phase[i], evall[i][j].fit[0].phase))

# Generate RMS vs Binsize plot
import correlated_noise as cn
rms     = []
stderr  = []
binsz   = []
binsz2  = []
normfactor = []
for i in range(numchan):
    foo     = cn.computeRMS(res[i], binstep=1)
    rms.append(foo[0])
    stderr.append(foo[1])
    binsz.append(foo[2])
    binsz2.append(binsz[i]*np.median(np.ediff1d(phase))*evall[0][0].period*86400)   # binsize in seconds
    normfactor.append(foo[1][0])

#normfactor = 1.
plt.figure(22, figsize=(6,4))
plt.clf()
plt.loglog(binsz2[0], rms[0]/normfactor[0], color='black', lw=1.5, label='Spectroscopic Light Curve RMS')    # our noise
for i in range(1,numchan):
    plt.loglog(binsz2[i], rms[i]/normfactor[i], color='black', lw=1.5)    # our noise
plt.loglog(binsz2[0], stderr[0]/normfactor[0], color='red', ls='-', lw=2, label='Standard Error') # expected noise
plt.xlim(binsz2[0][0], 5000)
plt.ylim(stderr[0][-1]/normfactor[0]/4., stderr[0][0]/normfactor[0]*2.)
plt.xlabel("Bin Size (seconds)", fontsize=12,labelpad=0)
plt.ylabel("Normalized RMS", fontsize=12,labelpad=0)
plt.xticks(size=12)
plt.yticks(size=12)
plt.legend(loc='upper right',fontsize=12)
plt.subplots_adjust(0.11,0.10,0.96,0.97)

plt.savefig('HD209458b-RMSvsBinSize-Spec.pdf')
plt.savefig('HD209458b-RMSvsBinSize-Spec.png')
"""


