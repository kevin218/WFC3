
import os, sys, copy
import manageevent as me
from run import *
from calc_scaleheight import *

# WASP 79b parameters
Rs = 1.64
mass = 0.9
radius = 1.7
temperature = 1959

filestr = []
filestr.append('/Users/kristin/Documents/STScI/KBS_WFC3-IR_Algorithm/2017-07-15-w1_spec_width_20/fitbghw_20/divideWLC_09072017/2017-09-18_08-26-test/')
filestr.append('/Users/kristin/Documents/STScI/KBS_WFC3-IR_Algorithm/2017-07-15-w1_spec_width_20/fitbghw_20/divideWLC_09272017/2017-09-28_16-34-test/')
filestr.append('/Users/kristin/Documents/STScI/KBS_WFC3-IR_Algorithm/2017-07-15-w1_spec_width_20/fitbghw_20/spectroscopicLC_08102017/2017-10-09_21-43-test/')
filestr.append('/Users/kristin/Documents/STScI/KBS_WFC3-IR_Algorithm/2017-07-15-w1_spec_width_20/fitbghw_20/spectroscopicLC_09272017/2017-10-09_21-31-test/')
#filestr.append('/Users/kristin/Documents/STScI/KBS_WFC3-IR_Algorithm/2017-07-15-w1_spec_width_20/fitbghw_20/spectroscopicLC_08102017/2017-10-03_21-15-test/')
#filestr.append('/Users/kristin/Documents/STScI/KBS_WFC3-IR_Algorithm/2017-07-15-w1_spec_width_20/fitbghw_20/spectroscopicLC_09272017/2017-09-28_16-28-test/')

plotColors = ['r','b','g','k']
saveFlag = 1

numspec = len(filestr)
wave_dw10 = []
waveerr_dw10 = []
depth_dw10 = []
deptherr_dw10 = []
wavelow_dw10 = []
wavehi_dw10 = []
linearLD_dw10 = []
transitH_dw10 = []

for j in range(len(filestr)):
    evTmp = w6Restore(filedir = filestr[j])
    numchan = len(evTmp)
    wave_tmp        = np.zeros(numchan)
    waveerr_tmp     = np.zeros(numchan)
    depth_tmp       = np.zeros(numchan)
    deptherr_tmp    = np.zeros(numchan)
    wavelow_tmp     = np.zeros(numchan)
    wavehi_tmp      = np.zeros(numchan)
    linearLD_tmp    = np.zeros(numchan)
    transitH_tmp    = np.zeros(numchan)

    for i in range(len(evTmp)):
        wave_tmp[i]     = (evTmp[i].wave_hi+evTmp[i].wave_low)/2.
        waveerr_tmp[i]  = (evTmp[i].wave_hi-evTmp[i].wave_low)/2.
        depth_tmp[i]    = 100*evTmp[i].fit[0].bestp[1]**2
        deptherr_tmp[i] = 200*evTmp[i].fit[0].bestp[1]*evTmp[i].fit[0].medianp[1,1]
        wavelow_tmp[i]  = evTmp[i].wave_low
        wavehi_tmp[i]   = evTmp[i].wave_hi
        linearLD_tmp[i] = evTmp[i].fit[0].bestp[6]
        transitH_tmp[i] = depth2scaleheight(depth_tmp[i], 1.13, Rs, mass, radius, temperature, mu=2.2)
    
    wave_dw10.append(wave_tmp)
    waveerr_dw10.append(waveerr_tmp)
    depth_dw10.append(depth_tmp)
    deptherr_dw10.append(deptherr_tmp)
    wavelow_dw10.append(wavelow_tmp)
    wavehi_dw10.append(wavehi_tmp)
    linearLD_dw10.append(linearLD_tmp)
    transitH_dw10.append(transitH_tmp)

#Plot spectrum
hanArr = []
fig = plt.figure()
plt.figure(1, figsize=(10,6))
plt.clf()
ax1 = fig.add_subplot(111)
ax2 = ax1.twinx()
plt.title('WASP-79b', size=16)
for j in range(numspec):
    hanTmp = ax1.errorbar(wave_dw10[j], depth_dw10[j], deptherr_dw10[j], waveerr_dw10[j], fmt=plotColors[j]+'o', ecolor=plotColors[j])
    ax2.plot(wave_dw10[j],transitH_dw10[j], plotColors[j]+'--')
    hanArr.append(hanTmp)
ax1.set_ylabel('Transit Depth (%)', size=14)
ax2.set_ylabel('Scale Height', size=14)
plt.xlabel('Wavelength ($\mu m$)', size=14)
plt.legend(hanArr, ['divideWLC10','divideWLC15','specLC10','specLC15'],loc='lower left')
plt.rc('grid',linestyle='--',color='gray')
ax1.grid()
#plt.ylim(0,200)
#plt.subplots_adjust(0.09,0.09,0.98,0.95)
if saveFlag:
    plt.savefig('WASP79b-WFC3-TrSpec_09302017.png')
    plt.savefig('WASP79b-WFC3-TrSpec_09302017.pdf')