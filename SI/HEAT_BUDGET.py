
import numpy as np
from scipy.stats import pearsonr

from grid_PAS import Grid as Grid_PAS

import matplotlib.pyplot as plt

import plotting as pt
import plotting_tools as ptt

import PAS_tools as ptools
import tools

import indices

from readData import *

#==========================================================

DATA_ROOT = 'data/'

HEAT_BUDGET_TIMESERIES = fff
if HEAT_BUDGET_TIMESERIES:

        print('To do on heat budget analysis: make isfw sum only contain heat fluxes going into the box area of integration. Could change box area of integration to more easily include ice shelves I average over. Wait for vvelth and uvelth from mitgcm run.')
        # Variables are defined, loaded and processed in block further below.   

        black = 'k'; blue = 'cornflowerblue'; red = 'indianred';
        green = 'seagreen'; grey = 'gray'

        nn = 48 

        year0 = 1979
        dt_month = 365.25 * 86400 / 12

        #==

        path = DATA_ROOT 

        # Load time
        t = np.load(path + 'PAS_time.npy')
        year = np.array(ptools.getDecimalTime(t))
        year = year - year[0] + year0
        year = ptools.windowAv(year, n=nn, av=False)

        #==

        # Empty list to add dictionaries to.
        data = []

        #==

        # Make edits/additions here.
        tot = 0
        totint = 0
        ##
        # LATERAL HEAT FLUX
        if True:
                tmp = np.load(path + 'heatBudget_uLateralFlux.npy')
                tmp += np.load(path + 'heatBudget_vLateralFlux.npy')

                # Subtract freezing temp contribution.
                #tmp -= -1.8*np.load(path + 'heatBudget_uCorrection.npy')
                #tmp -= -1.8*np.load(path + 'heatBudget_vCorrection.npy')
                tmp = ptools.detrend(tmp)
                tmp = ptools.windowAv(tmp, n=nn)

                tmpint = dt_month * np.cumsum(tmp)
                tmpint = ptools.detrend(tmpint)
                tot += tmp
                totint += tmpint
                label = r'Lateral heat flux'
                dct = {'data':tmp, 'dataint':tmpint, 'time':year, 'timeint':year, 'label':label, 'color':red}
                data.append(dct)
        ##

        ##
        # ICE-SHELF FW HEAT RELEASE
        if True:
                tmp = np.load(path + 'heatBudget_isfw.npy')
                tmp = ptools.detrend(tmp)
                tmp = ptools.windowAv(tmp, n=nn)
                tmpint = dt_month * np.cumsum(tmp)
                tmpint = ptools.detrend(tmpint)
                tot += tmp
                totint += tmpint
                label = r'Ice-shelf melt latent heat release'
                dct = {'data':tmp, 'dataint':tmpint, 'time':year, 'timeint':year, 'label':label, 'color':green}
                data.append(dct)
        ##

        ##
        # SURFACE HEAT FLUX
        if True:
                tmp = np.load(path + 'heatBudget_oceQnet.npy')
                tmp = ptools.detrend(tmp)
                tmp = ptools.windowAv(tmp, n=nn)
                tmpint = dt_month * np.cumsum(tmp)
                tmpint = ptools.detrend(tmpint)
                tot += tmp
                totint += tmpint
                label = r'Surface heat flux'
                dct = {'data':tmp, 'dataint':tmpint, 'time':year, 'timeint':year, 'label':label, 'color':blue}
                data.append(dct)
        ##
        ##
        # DEPTH-INTEGRATED TEMPERATURE
        if True:
                tmp = np.load(path + 'heatBudget_heatContent.npy')
                tmp = ptools.detrend(tmp)
                tmpint = ptools.windowAv(tmp, n=nn)
                tmp = (tmp[1:] - tmp[:-1]) / dt_month
                tmp = ptools.windowAv(tmp, n=nn)
                tmpint = ptools.detrend(tmpint)
                #tmpint = dt_month * np.cumsum(tmp) 
                #tot += tmp
                #totint += tmpint
                label = r'Heat content change'
                dct = {'data':tmp, 'dataint':tmpint, 'time':year[1:], 'timeint':year, 'label':label, 'color':'black'}
                data.append(dct)
                ERROR = np.zeros(year.shape[0]); ERROR[1:] = tot[1:] - tmp; ERROR[0] = ERROR[1]
                dct = {'data':ERROR, 'dataint':totint-tmpint, 'time':year, 'timeint':year, 'label':'Error', 'color':'grey'}
                data.append(dct)
        ##

        #==

        fs = 14

        # PLOT
        plt.figure(figsize=(15,6))
        for var in data:
                data_tmp = var['data']
                label = var['label']
                #aaa
                if label == r'Lateral heat flux':
                        data_tmp += ERROR / 2
                elif label == 'Error':
                        data_tmp -= ERROR / 2
                #bbb
                plt.plot(var['time'], 1.e-12*data_tmp, label=label, color=var['color'])
                ucstr = ''
        #plt.plot(year, data[0]['data']+data[1]['data'], color='orange', label='Lat flux - ISFW')
        plt.xlim(np.min(var['time']), np.max(var['time']))
        plt.grid()
        plt.title(r'On-shelf heat budget timeseries ($10^{12}$ W)', fontsize=fs+2)
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)
        plt.legend(prop={'size':fs-1})
        plt.tight_layout()
        plt.savefig('heat_budget_timeseries')
        plt.show()


