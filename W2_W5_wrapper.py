%matplotlib inline

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import multiprocessing as mp
import manageevent as me
import time, optspex, sys
from astropy.io import fits
import imp

# SPectrum Analysis, Reduction, and Calibration (SPARC) pipeline
# Written by Kevin B. Stevenson in June of 2012
import w1reduce as w1
import w2reduce as w2
import w3lc as w3
import w4ld as w4
import w5lc as w5

#eventname = 'hat38v1'
eventname = 'wasp79b'
eventdir   = '2017-07-09-w1'
# eventdir   = '2017-05-17-w2-256'

import imp

madVarList = ['p7thresh']
varColors = ['m']
#madVarList = ['diffthresh','p3thresh','p5thresh','p7thresh','specwidth','fitbghw']
#varColors = ['b','g','c','m','r','k']
madVarRangeList = [[3,5,7,10,12,15,20]]
#madVarRangeList = [[5,7,10,12,15],[3,4,5,6,7],[5,7,10,12,15],[5,7,10,12,15],[10,14,17,20,25],[10,14,17,20,25]]
isplots = 1

for j in range(len(madVarList)):
    madVariable = madVarList[j]
    madVarRange = madVarRangeList[j]
    for i in range(len(madVarRange)):
        madVarSet = madVarRange[i]
        print(madVariable + ' = ' + str(madVarSet))
        # import w2reduce as w2
        eventdir   = '2017-07-15-w1_'+madVariable+'_'+str(madVarSet)
        ev=w2.reduceWFC3(eventname, eventdir, madVariable, madVarSet, isplots=1)
        # w5 requires w2 output
        w5.lcWFC3(eventname, eventdir, 1, madVariable, madVarSet, isplots=1, expand=1)
    
    f2 = open('W5_MAD_'+ madVariable +'_1D.txt','r')
    madvar2,mad2_flux0,mad2_flux1 = np.loadtxt('W5_MAD_'+ madVariable +'_1D.txt',unpack=True,delimiter=',')
    #lines2 = f2.readlines()

    f5 = open('W5_MAD_'+ madVariable +'.txt','r')
    madvar5,mad5 = np.loadtxt('W5_MAD_'+ madVariable +'.txt',unpack=True,delimiter=',')

    if isplots:
        plt.figure(3000, figsize=(10,8))
        plt.clf()
        plt.subplot(211)
        plt.plot(madvar2, mad2_flux0, "o", color=varColors[j])
        plt.title("W5 1D MAD vs "+madVariable)
        plt.ylabel('MAD')
        plt.xlabel(madVariable)
        plt.grid()
        plt.tight_layout()

        plt.subplot(212)
        plt.plot(madvar5, mad5, "o", color=varColors[j])
        plt.title("W5 MAD vs "+madVariable)
        plt.ylabel('MAD')
        plt.xlabel(madVariable)
        plt.grid()
        plt.tight_layout()
        plt.savefig(eventdir+'/figs/MAD_vs_'+madVariable+'.png')