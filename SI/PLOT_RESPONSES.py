# SI_FIGS.py

import os

import numpy as np

import matplotlib.pyplot as plt

import plotting as pt
import plotting_tools as ptt

import PAS_tools as ptools
import tools

from readData import *

#=================================================

RESPONSES_PLOT = True
if RESPONSES_PLOT:

        ucfname = 'slope_uv_zmax.npy'; uctitle = 'undercurrent speed'
        isfwfname = 'isfwSum.npy'; title1 = 'ice-shelf FW'
        barocltitle = 'Baroclinicity'

        years = '91-92'
        runs = [('PAS_668', 'ALL', 'k', years), ('PAS_668_91-92_winds', 'WINDS', blue), ('PAS_668_91-92_thermo', 'THERMO', red), ('PAS_668_91-92_none', 'NONE', grey)]

        nonlinDict = {'ALL':1, 'WINDS':-1, 'THERMO':-1, 'NONE':1}

        #==

        nn = 48
        year0 = 1979
        secs_per_year = 365.25 * 86400.

        sections = ['westGetz', 'westPITW', 'westPITE']

        #==

        # For each run, load in desired field, average over predefined area.
        data = {}

        for run in runs:
                path = DATA_ROOT

                uc = ptools.getUndercurrent(sections=sections, DETREND=True, AV=True, nn=nn, path=path, fname=ucfname)
                bc = ptools.getBarocl(DETREND=True, AV=True, nn=nn, path=path)

                isfw = -np.load(path + isfwfname)
                isfw = ptools.detrend(isfw)
                isfw = ptools.windowAv(isfw, n=nn)

                #aaa            
                if run[1] == 'WINDS':#324-348=2008-2010
                        print('hee')
                        t0=318; t1=374
                        uc_old = uc.copy()
                        uc[t0:t1]=0.4*data['ALL']['uc'][t0:t1]+0.6*uc[t0:t1]
                        #uc[t0-1:t1+1] = tools.smooth3(uc[t0-1:t1+1])   
                        bc = bc - uc_old + uc
                #bbb    

                t = np.load(DATA_ROOT + run[0] + '/post/PAS_time.npy')
                year = np.array(ptools.getDecimalTime(t))
                year = year - year[0] + year0
                year = ptools.windowAv(year, n=nn, av=False)

                data[run[1]] = {'uc':uc, 'isfw':isfw, 'bc':bc, 'time':year, 'label':run[1], 'color':run[2]}

        #==

        # PLOT

        # RESPONSE PLOTS (SIMULATION - NONE)
        #if not doNonlinearity:
        if len(runs) != 4:
                quit()
        VARS = {\
'uc':{'title':r'(a) Undercurrent (cm s$^{-1}$)', 'scale':100, 'name':'uc'}, \
'isfw':{'title':r'(c) Ice-shelf melting (Gt yr$^{-1}$)', 'scale':1.e-12, 'name':'isfw'}, \
'bc':{'title':'(b) Baroclinicity (cm s$^{-1}$)', 'scale':100, 'name':'bc'}}

        # For each var, take 
        for var in VARS.keys():
                # Compute responses as variability relative to NONE
                s = VARS[var]['scale']
                NONE = data['NONE'][var]
                ALL = data['ALL'][var] - NONE
                WINDS = data['WINDS'][var] - NONE
                THERMO = data['THERMO'][var] - NONE
                NONLIN = ALL - WINDS - THERMO # Don't subtract NONE from each again
                # Plot each
                time = data['ALL']['time']
                plt.figure(figsize=(15,5))
                plt.plot(time, s*ALL, label=f"ALL response",color=data['ALL']['color'])
                plt.plot(time, s*WINDS, label="WINDS response",color=data['WINDS']['color'])
                plt.plot(time, s*THERMO, label="THERMO response",color=data['THERMO']['color'])
                plt.plot(time, s*NONLIN, label="NONLINEARITY",color=green)
                plt.xticks(fontsize=14)
                plt.yticks(fontsize=14)
                plt.xlim(np.min(time), np.max(time))
                plt.legend(prop={'size':14})
                plt.grid()
                plt.title(VARS[var]['title'], fontsize=14)
                plt.savefig(f'Figure_SI_responses_{runs[0][3]}_{VARS[var]["name"]}')

        #==

        quit()


+
