#! /usr/bin/env python

# $Author: kevin $
# $Revision: 540 $
# $Date: 2011-08-13 00:08:23 -0400 (Sat, 13 Aug 2011) $
# $HeadURL: file:///home/esp01/svn/code/python/pipeline/trunk/run.py $
# $Id: run.py 540 2011-08-13 04:08:23Z kevin $

"""
#TO EXECUTE DIRECTLY FROM BASH, USE THE FOLLOWING SYNTAX:
./run.py p5idl filedir
./run.py p5 filedir topdir clip
./run.py w6 newdir filedir topdir clip idl
./run.py w7 filedir topdir idl islak
./run.py p8 filedir topdir idl eclphase
./run.py p9 filedir topdir idl
#Arguments after p# are optional but must remain in order.
#   eg. To specify topdir with w7, you MUST include filedir;
#       however, specifying filedir does NOT require topdir or idl.
#newdir:    String to append to the default directory name.
            A new directory is ALWAYS created in BASH mode.
#filedir:   Location of the savefile to restore, default is '..' for w6 and '.' for others
#topdir:    Specify POET top directory if it's not a parent of cwd.
#           To specify default directory, type None.
#clip:      Clip data set to model and plot only a portion of the entire light curve.
#           Good for modeling a transit or eclipse in an around-the-orbit data set.
#           Syntax is start:end (eg. 0:70000 or -60500:None), defult is None for full data set.
#idl:       False (default), set True when using IDL photometry.
#islak:     True (default), set False to not load 'allknots' into memory and skip Figs 701 & 702.
#Figures will NOT display in BASH mode, but they will still save.
"""

from __future__ import print_function
import os, sys, copy
os.environ['OMP_NUM_THREADS']='6'
sys.path.append('../')
from ..lib import manageevent as me
from run import *
ancildir      = 'ancil/'
modeldir      = 'ancil/modelparams/'
isinteractive = True       #Set False and don't use '-pylab' to turn off displaying plots in interactive mode.

def ancil(ev):
    #ev = [ev[0],copy.deepcopy(ev[0])]
    #ev[0].eventname = 'wa063bph1_0'
    #ev[1].eventname = 'wa063bph1_1'
    for i in range(len(ev)):
        #j       = 0
        j       = i % 2
        iscan   = np.where(ev[i].scandir == j)[0]
        #print('scandir = '+str(ev[i].scandir))
        ev[i].aplev = np.sum(ev[i].photflux[0,iscan],axis=1)
        ev[i].aperr = np.sqrt(np.sum(ev[i].photfluxerr[0,iscan]**2,axis=1))
        ev[i].good  = np.ones(len(ev[i].aplev))
        ev[i].good[np.where(ev[i].orbitnum[iscan] == 0)[0]] = 0 #Remove first orbit
        #ev[i].good[:12] = 0
        ev[i].bjdtdb= ev[i].bjdtdb[iscan]
        ev[i].framenum= ev[i].framenum[iscan]
        ev[i].batchnum= ev[i].batchnum[iscan]
        ev[i].orbitnum= ev[i].orbitnum[iscan]
        ev[i].phase = ((ev[i].bjdtdb-ev[i].ephtime) % ev[i].period)/ev[i].period
        ev[i].npos  = 1
        ev[i].pos   = np.zeros(ev[i].n_files)
        ev[i].paramsfile   = ancildir + ev[i].eventname + 'params.py'
        ev[i].initvalsfile = ancildir + ev[i].eventname + '-initvals.txt'
        # Copy eg00params
        if os.path.isfile(ev[i].paramsfile) == False:
            shutil.copy(ancildir + 'eg00params.py', ev[i].paramsfile)
        # Copy eg00-initvals
        if os.path.isfile(ev[i].initvalsfile) == False:
            shutil.copy(ancildir + 'eg00-initvals.txt', ev[i].initvalsfile)
        if j == 0:
            ev[i].good[np.where(ev[i].framenum == 0)[0]] = 0
    #ev[0].good[[46,47]] = 0
    #ev[0].good[[27]] = 0
    return ev

#USE THE COMMANDS BELOW WHEN IN INTERACTIVE MODE
#ONLY COPY AND PASTE THE LINES THAT YOU NEED
def interactive():
    #Restore event

    #ev              = me.loadevent('d-wa012-p3')
    #ev[0].eventname = 'wa012bg1'
    #ev = w4Restore(filedir='..', fname='d-hd209bhs1_1125_1650-w3-WHITE.dat')
    ev = w4Restore(filedir='..', fname='d-wa063bph1_1125_1650-w3-WHITE.dat')
    #ev = w4Restore(filedir='..', fname=None)
    #ev = None
    ev = ancil(ev)

    #Run w6model
    #ev:        If not specified, event will be restored using poetRestore or p5idlRestore.
    #newdir:    Creates a new model directory when True (default).
    #           Can set to False or a string that appends to the default directory name.
    #filedir:   If event is not specified, location of the savefile to restore.
    #idl:       False (default), set True when using IDL photometry.
    w6model(ev)
    w6model(ev, newdir=False)
    w6model(ev, newdir='test')

    #Run w7anal, p8tables, and p9figs
    #event:     If not specified, event will be restored using w6Restore.
    #filedir:   If event is not specified, location of the savefile to restore.
    #topdir:    Specify POET top directory if it's not a parent of cwd.
    #idl:       False (default), set True when using IDL photometry.
    #islak:     True (default), set False to not load 'allknots' into memory and skip Figs 701 & 702 (w7anal only).
    #eclphase:  Nominal eclipse phase, default is 0.5 (p8tables only).
    w7anal(ev)
    p8tables(ev)
    p9figs(ev)

    #Restore after w6model or w7anal, if necessary
    #You can use os.listdir('.') to find the desired model directory.
    #filedir:   Location of the savefile to restore.
    #topdir:    Specify POET top directory if it's not a parent of cwd.
    #idl:       False (default), set True when using IDL photometry.
    ev = w6Restore(filedir = '2016-10-14_17:09-saveOutput')
    ev = w6Restore(filedir='.')
    ev = w7Restore(filedir='.')

    #np.savez('WASP63b-WhiteLC.npz', bjdtdb=ev[0].bjdtdb, flux=ev[0].aplev, err=ev[0].aperr)

    return

########################################################
#                                                      #
#   DO NOT RUN CODE INTERACTIVELY BEYOND THIS POINT!   #
#                                                      #
########################################################

def w4Restore(filedir='..', fname=None, topdir=None, clip=None):
    import shutil
    global numevents

    files    = []
    event    = []
    filename = ''
    if fname == None:
        for fname in os.listdir(filedir):
            if (fname.endswith("-w3.dat")):
                files.append(fname[:-4])
    else:
        files.append(fname.rstrip('.dat'))
    files.sort()
    numevents = len(files)
    if numevents == 0:
        print('Cannot find any files to restore.')
        #event = ancil(None)
        return event
    for i in np.arange(numevents):
        #Load event
        event.append(me.loadevent(filedir+'/'+files[i]))
        print('Finished loading: ' + event[i].eventname)
        filename = ''.join((filename,event[i].eventname))
        event[i].ancildir   = ancildir
    #On initial setup, rename eg00params and eg00-initvals
    '''
    for i in np.arange(numevents):
        event[i].paramsfile   = ancildir + event[i].eventname + 'params.py'
        event[i].initvalsfile = ancildir + event[i].eventname + '-initvals.txt'
        # Copy eg00params
        if os.path.isfile(event[i].paramsfile) == False:
            shutil.copy(ancildir + 'eg00params.py', event[i].paramsfile)
        # Copy eg00-initvals
        if os.path.isfile(event[i].initvalsfile) == False:
            shutil.copy(ancildir + 'eg00-initvals.txt', event[i].initvalsfile)
    '''
    return event

#RUN w6model
def w6model(ev=None, newdir=True, filedir='..', topdir=None, clip=None, idl=False):
    import datetime
    import w6model as w6
    #reload(w6)
    global nummodels, isinteractive
    if ev == None or len(ev) == 0:
        print("Event object is empty.")
        ev = w4Restore(filedir='..', fname=None)
        ev = ancil(ev)


    if newdir == True:
        #CREATE NEW DIRECTORY FOR MODEL RUN
        modeldir = datetime.datetime.strftime(datetime.datetime.today(), "%Y-%m-%d_%H:%M")
        os.mkdir(modeldir)
    elif type(newdir) == str:
        #CREATE NEW DIRECTORY FOR MODEL RUN WITH CUSTOM NAME
        modeldir = datetime.datetime.strftime(datetime.datetime.today(), "%Y-%m-%d_%H:%M")+"-"+newdir
        os.mkdir(modeldir)
    else:
        try:
            modeldir = ev[0].modeldir
        except:
            print("Model directory has not been specified.")
            return

    #RELOAD MODEL PARAMETERS
    numevents = len(ev)
    nummodels = np.zeros(numevents,dtype=int)
    for j in range(numevents):
        exec('import ' + ev[j].eventname + 'params as op' + str(j))
        #exec("reload(op" + str(j) + ")")
        ev[j].params = readeventhdf.initParams()
        exec("op" + str(j) + ".modelparams(ev["+ str(j) + "].params)")
        #exec("ev[" + str(j) + "].params = op" + str(j) + ".params")
        ev[j].fit = []
        ev[j].modeldir = modeldir
        nummodels[j] = len(ev[j].params.model)
        if j > 0 and nummodels[j] != nummodels[j-1]:
            print("WARNING: Number of models in each event does not match.")

    #INITIALIZE OUTPUT TYPE: stdout, A FILE OBJECT OR A FILE
    printout = printoutput.init(ev[0].params.printout, ev)

    #Execute rundmc
    for i in range(nummodels.min()):
        w6.rundmc(ev, i, printout, isinteractive)
        if ev[0].params.savedata:
            for j in range(numevents):
                w6Save(ev[j])

    #PRINT PARAMETERS USED FOR COMPARISON
    print("\nBest-fit eclipse depths or transit radius ratios with errors:", file=printout)
    for j in range(numevents):
        ev[j].minbic = np.inf     #Minimum BIC value of all fits for one event.
        print(ev[j].eventname, file=printout)
        for i in range(len(ev[j].fit)):
            if hasattr(ev[j].fit[i].i,'depth'):
                print(ev[j].fit[i].bestp  [ev[j].fit[i].i.depth  ],
                      ev[j].fit[i].medianp[ev[j].fit[i].i.depth,1],
                      ev[j].fit[i].saveext, file=printout)
            if hasattr(ev[j].fit[i].i,'depth2'):
                print(ev[j].fit[i].bestp  [ev[j].fit[i].i.depth2  ],
                      ev[j].fit[i].medianp[ev[j].fit[i].i.depth2,1],
                      ev[j].fit[i].saveext, file=printout)
            if hasattr(ev[j].fit[i].i,'depth3'):
                print(ev[j].fit[i].bestp  [ev[j].fit[i].i.depth3  ],
                      ev[j].fit[i].medianp[ev[j].fit[i].i.depth3,1],
                      ev[j].fit[i].saveext, file=printout)
            if hasattr(ev[j].fit[i].i,'trrprs'):
                print(ev[j].fit[i].bestp  [ev[j].fit[i].i.trrprs  ],
                      ev[j].fit[i].medianp[ev[j].fit[i].i.trrprs,1],
                      ev[j].fit[i].saveext, file=printout)
            if hasattr(ev[j].fit[i].i,'trrprs2'):
                print(ev[j].fit[i].bestp  [ev[j].fit[i].i.trrprs2  ],
                      ev[j].fit[i].medianp[ev[j].fit[i].i.trrprs2,1],
                      ev[j].fit[i].saveext, file=printout)
            if hasattr(ev[j].fit[i].i,'trqrprs'):
                print(ev[j].fit[i].bestp  [ev[j].fit[i].i.trqrprs  ],
                      ev[j].fit[i].medianp[ev[j].fit[i].i.trqrprs,1],
                      ev[j].fit[i].saveext, file=printout)
            if hasattr(ev[j].fit[i].i,'rprs'):
                print(ev[j].fit[i].bestp  [ev[j].fit[i].i.rprs  ],
                      ev[j].fit[i].medianp[ev[j].fit[i].i.rprs,1],
                      ev[j].fit[i].saveext, file=printout)
            if hasattr(ev[j].fit[i].i,'rprs2'):
                print(ev[j].fit[i].bestp  [ev[j].fit[i].i.rprs2  ],
                      ev[j].fit[i].medianp[ev[j].fit[i].i.rprs2,1],
                      ev[j].fit[i].saveext, file=printout)
            ev[j].minbic = np.min((ev[j].minbic,ev[j].fit[i].bic))

    print("\n     S/N      SDNR     \xce\x94BIC       MODEL   NUMIT  BIN_SZ(y,x)   MinNumPts", file=printout)
    #Delta = '\xce\x94' in utf-8
    for j in range(numevents):
        print(ev[j].eventname, file=printout)
        minbic = ev[j].minbic
        for i in range(len(ev[j].fit)):
            try:
                sdnr  = ev[j].fit[i].sdnr
                bic   = ev[j].fit[i].bic - minbic
                model = ev[j].fit[i].saveext
                numit = ev[j].params.numit[1]
                if len(ev[j].params.ystep) == len(ev[j].params.model):
                    ystep, xstep = ev[j].params.ystep[i], ev[j].params.xstep[i]
                else:
                    ystep, xstep = ev[j].params.ystep[0], ev[j].params.xstep[0]
                if len(ev[j].params.minnumpts) == len(ev[j].params.model):
                    minnumpts = ev[j].params.minnumpts[i]
                else:
                    minnumpts = ev[j].params.minnumpts[0]
                if hasattr(ev[j].fit[i].i,'depth'):
                    snr   = ev[j].fit[i].bestp  [ev[j].fit[i].i.depth  ] / \
                            ev[j].fit[i].medianp[ev[j].fit[i].i.depth,1]
                    print('%8.4f %9.7f %8.1f %11s %7.1e %6.3f,%5.3f %4.0f' %
                         (snr, sdnr, bic, model, numit, ystep, xstep, minnumpts), file=printout)
                if hasattr(ev[j].fit[i].i,'depth2'):
                    snr   = ev[j].fit[i].bestp  [ev[j].fit[i].i.depth2  ] / \
                            ev[j].fit[i].medianp[ev[j].fit[i].i.depth2,1]
                    print('%8.4f %9.7f %8.1f %11s %7.1e %6.3f,%5.3f %4.0f' %
                         (snr, sdnr, bic, model, numit, ystep, xstep, minnumpts), file=printout)
                if hasattr(ev[j].fit[i].i,'depth3'):
                    snr   = ev[j].fit[i].bestp  [ev[j].fit[i].i.depth3  ] / \
                            ev[j].fit[i].medianp[ev[j].fit[i].i.depth3,1]
                    print('%8.4f %9.7f %8.1f %11s %7.1e %6.3f,%5.3f %4.0f' %
                         (snr, sdnr, bic, model, numit, ystep, xstep, minnumpts), file=printout)
                if hasattr(ev[j].fit[i].i,'trrprs'):
                    snr   = ev[j].fit[i].bestp  [ev[j].fit[i].i.trrprs  ] / \
                            ev[j].fit[i].medianp[ev[j].fit[i].i.trrprs,1]
                    print('%8.4f %9.7f %8.1f %11s %7.1e %6.3f,%5.3f %4.0f' %
                         (snr, sdnr, bic, model, numit, ystep, xstep, minnumpts), file=printout)
                if hasattr(ev[j].fit[i].i,'trrprs2'):
                    snr   = ev[j].fit[i].bestp  [ev[j].fit[i].i.trrprs2  ] / \
                            ev[j].fit[i].medianp[ev[j].fit[i].i.trrprs2,1]
                    print('%8.4f %9.7f %8.1f %11s %7.1e %6.3f,%5.3f %4.0f' %
                         (snr, sdnr, bic, model, numit, ystep, xstep, minnumpts), file=printout)
                if hasattr(ev[j].fit[i].i,'trqrprs'):
                    snr   = ev[j].fit[i].bestp  [ev[j].fit[i].i.trqrprs  ] / \
                            ev[j].fit[i].medianp[ev[j].fit[i].i.trqrprs,1]
                    print('%8.4f %9.7f %8.1f %11s %7.1e %6.3f,%5.3f %4.0f' %
                         (snr, sdnr, bic, model, numit, ystep, xstep, minnumpts), file=printout)
                if hasattr(ev[j].fit[i].i,'rprs'):
                    snr   = ev[j].fit[i].bestp  [ev[j].fit[i].i.rprs  ] / \
                            ev[j].fit[i].medianp[ev[j].fit[i].i.rprs,1]
                    print('%8.4f %9.7f %8.1f %11s %7.1e %6.3f,%5.3f %4.0f' %
                         (snr, sdnr, bic, model, numit, ystep, xstep, minnumpts), file=printout)
                if hasattr(ev[j].fit[i].i,'rprs2'):
                    snr   = ev[j].fit[i].bestp  [ev[j].fit[i].i.rprs2  ] / \
                            ev[j].fit[i].medianp[ev[j].fit[i].i.rprs2,1]
                    print('%8.4f %9.7f %8.1f %11s %7.1e %6.3f,%5.3f %4.0f' %
                         (snr, sdnr, bic, model, numit, ystep, xstep, minnumpts), file=printout)
            except:
                print("Error calculating values. %13s" % ev[j].fit[i].saveext, file=printout)
    printoutput.close(printout)
    if isinteractive == False:
        plt.close('all')
    else:
        plt.show()
    return

#SAVE event AFTER w6model
def w6Save(ev):
    cwd = os.getcwd().split("/")
    if cwd[-1] == ev.modeldir:
        savefile  = "d-" + ev.eventname + "-6model.dat"
    else:
        savefile  = ev.modeldir + "/d-" + ev.eventname + "-6model.dat"
    '''
    handle    = open(savefile,'wb')
    pickle.dump(ev, savefile)
    handle.close()
    pickle.dump(ev, handle)
    '''
    with open(savefile,'wb') as pickle_file:
        pickle.dump(ev, pickle_file)
    return

#RESTORE SAVEFILE AFTER w6model
def w6Restore(filedir='.', topdir=None, idl=False):
    from ..lib import np_unpickle as unpic #added
    global numevents, nummodels
    '''
    if idl == False:
        #Append system path
        if topdir == None:
            r = os.getcwd().split("/")
            topdir = "/".join(r[:r.index("run")])
        sys.path.append(topdir + '/lib/')
    '''
    loadfile = []
    for fname in os.listdir(filedir):
        if (fname.endswith("6model.dat")):
            loadfile.append(filedir+"/"+fname)
    loadfile.sort()
    #print("Loading files:", loadfile)
    #loadfile = [loadfile[-1]]   #***EDIT THIS LINE MANUALLY***
    numevents = len(loadfile)
    if numevents == 0:
        print('Cannot find any files to restore.')
        return []
    ev    = []
    for lfile in loadfile:
        print("Loading " + lfile)
        handle      = open(lfile, 'rb')
        try:
            ev.append(pickle.load(handle))
            #ev.append(pickle.load(lfile))
        except:
            ev.append(unpic.unpickle_old_pyfits(lfile))
        #handle.close()
    nummodels = np.zeros(numevents,dtype=int)
    for j in range(numevents):
        nummodels[j] = len(ev[j].params.model)
    return ev

#RUN 7anal
def w7anal(ev=None, filedir='.', topdir=None, idl=False, islak=True):
    import w7anal as w7
    
    global numevents, nummodels, isinteractive
    if ev == None:
        ev = w6Restore(filedir, topdir, idl)
    cwd = os.getcwd().split("/")
    if cwd[-1] == ev[0].modeldir:
        os.chdir('..')
    printout = printoutput.init(ev[0].params.printout, ev)
    for j in range(numevents):
        print("\n" + ev[j].eventname, file=printout)
        for i in range(nummodels.min()):
            print("\nCurrent model = " + str(ev[j].params.model[i]), file=printout)
            w7.stdanal(ev[j], ev[j].fit[i], j*nummodels.min()+i, printout, islak)
            w7Save(ev[j])
    printoutput.close(printout)
    if isinteractive == False:
        plt.close('all')
    else:
        plt.show()
    return

#SAVE event AFTER w7anal
def w7Save(ev):
    cwd = os.getcwd().split("/")
    if cwd[-1] == ev.modeldir:
        savefile  = "d-" + ev.eventname + "-7anal.dat"
    else:
        savefile  = ev.modeldir + "/d-" + ev.eventname + "-7anal.dat"
    handle    = open(savefile, 'w')
    pickle.dump(ev, handle)
    handle.close()
    return

#RESTORE SAVEFILE AFTER w7anal
def w7Restore(filedir='.', topdir=None, idl=False):
    global numevents, nummodels
    if idl == False:
        #Append system path
        if topdir == None:
            r = os.getcwd().split("/")
            topdir = "/".join(r[:r.index("run")])
        sys.path.append(topdir + '/lib/')
    loadfile = []
    for fname in os.listdir(filedir):
        if (fname.endswith("7anal.dat")):
            loadfile.append(filedir+"/"+fname)
    loadfile.sort()
    numevents = len(loadfile)
    if numevents == 0:
        print('Cannot find any files to restore.')
        return []
    ev    = []
    for lfile in loadfile:
        print("Loading " + lfile)
        handle      = open(lfile, 'r')
        try:
            ev.append(pickle.load(handle))
        except:
            ev.append(unpic.unpickle_old_pyfits(lfile))
        handle.close()
    nummodels = np.zeros(numevents,dtype=int)
    for j in range(numevents):
        nummodels[j] = len(ev[j].params.model)
    return ev

#RUN p8tables
def p8tables(ev=None, filedir='.', topdir=None, idl=False, eclphase=0.5):
    from ..lib import p8tables as p8
    
    global numevents, nummodels
    if ev == None:
        ev = w7Restore(filedir, topdir, idl)
    cwd = os.getcwd().split("/")
    if cwd[-1] == ev[0].modeldir:
        os.chdir('..')
    printout = printoutput.init(ev[0].params.printout, ev)
    for j in range(numevents):
        print("\n" + ev[j].eventname, file=printout)
        ev[j].meanphase = eclphase
        for i in range(nummodels.min()):
            print("\nCurrent model = " + str(ev[j].params.model[i]), file=printout)
            p8.tables(ev[j], i, printout)
    printoutput.close(printout)
    return

#RUN p9figs
def p9figs(ev=None, filedir='.', topdir=None, idl=False):
    from ..lib import p9figs as p9
    
    global numevents, nummodels, isinteractive
    if ev == None:
        ev = w7Restore(filedir, topdir, idl)
    cwd = os.getcwd().split("/")
    if cwd[-1] == ev[0].modeldir:
        os.chdir('..')
    printout = printoutput.init(ev[0].params.printout, ev)
    for j in range(numevents):
        print("\n" + ev[j].eventname, file=printout)
        for i in range(nummodels.min()):
            print("\nCurrent model = " + str(ev[j].params.model[i]), file=printout)
            p9.figs  (ev[j], i, j*nummodels.min()+i)
    printoutput.close(printout)
    if isinteractive == False:
        plt.close('all')
    else:
        plt.show()
    return

import numpy as np
import matplotlib.pyplot as plt
import _pickle as pickle
#import cPickle as pickle
import shutil
from ..lib import printoutput, readeventhdf
from ..lib import np_unpickle as unpic
sys.path.append(ancildir)
sys.path.append('../'+ancildir)

#MAIN IS CALLED WHEN EXECUTING DIRECTLY FROM BASH
def main(args):
    length = len(args)
    #Run w6model
    if args[1] == 'w6':
        w6model(None, *args[2:])
    #Run w7anal
    elif args[1] == 'w7':
        w7anal(None, *args[2:])
    #Run p8tables
    elif args[1] == 'w8':
        p8tables(None, *args[2:])
    #Run p9figs
    elif args[1] == 'w9':
        p9figs(None, *args[2:])
    else:
        print("Unrecognized function.")
    return 0

#CALLS main IN BASH MODE THEN EXITS CLEANLY
if __name__ == '__main__':
    isinteractive = False
    sys.exit(main(sys.argv))
