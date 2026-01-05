# SI_FIGS.py

import os

import numpy as np

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

black = 'k'; blue = 'cornflowerblue'; red = 'indianred'; green = 'seagreen'; grey = 'gray'
cmapMean = 'plasma'; cmapAnom = 'coolwarm'; white = 'w'; cyan = 'cyan'; darkgrey = 'dimgrey'

# Time constants
secs_per_year = 86400. * 365.25
secs_per_month = secs_per_year / 12.

# Factor for converting kg/m2/s to m/yr.
fwConversion = 1.e-3 * secs_per_year

# Lats for Amundsen Sea subregion.
EASlons = [235, 262] 
EASlats = [-75.5, -70]

LATS = [-75.5, -68]
LONS = [225, 262]

# Area of integration. Set longitudinal ranges on/off the shelf.
#lonsN=[225., 245]; lonsS=[233,258.5]; Nrange=58
lonsN=[230., 245]; lonsS=[233,258.5]; Nrange=58
		
#==========================================================

ttt = True
fff = False

#==========================================================

PLOT_NONES = fff
if PLOT_NONES:

	ucfname = 'slope_uv_zmax.npy'; uctitle = 'undercurrent speed'
	isfwfname = 'isfwSum.npy'; isfwtitle = 'ice-shelf FW'

	refrun = 'PAS_668'
	
	YEARS = ['79-80', '83-84', '87-88', '91-92', '95-96', '99-00', '03-04', '07-08', '11-12', '15-16', '19-20']
	runs = [f'PAS_668_{year}_none' for year in YEARS]
	runs.append(refrun)

	# Factor for converting area-integrated FW fluxes into Gt/yr
	s = 1.e-12 * 365.25 * 86400 

	#==
	
	nn = 48

	year0 = 1979
	
	sections = ['westGetz', 'westPITW', 'westPITE']
	
	#==

	# For each run, load in desired field (e.g., FW or SIheff), average over predefined area.
	data = {}
	
	# Do reference run first
	
	ucmax = 0; isfwmax = 0
	for run in runs:
		path = DATA_ROOT 
		
		uc = ptools.getUndercurrent(sections=sections, DETREND=False, AV=True, nn=nn, path=path, fname=ucfname)
		ucmax = max(ucmax, np.max(np.abs(uc)))
		
		isfw = -1.e-12*np.load(path + isfwfname)
		isfw = ptools.windowAv(isfw, n=nn)
		isfwmax = max(isfwmax, np.max(np.abs(isfw)))
	
		t = np.load(path+'PAS_time.npy')
		year = np.array(ptools.getDecimalTime(t))
		year = year - year[0] + year0
		year = ptools.windowAv(year, n=nn, av=False)
		
		if run == refrun:
			label = 'ALL'
		else:
			label = str(int(1e2*(20-np.floor(int(run[8:10])/50))+int(run[8:10])))+'-'+str(int(1e2*(20-np.floor(int(run[11:13])/50))+int(run[11:13])))

		data[run] = {'uc':uc, 'isfw':isfw, 'time':year, 'label':label}
	
	#==
	
	# PLOT
	
	# 0. Undercurrent
	nonlin = 0
	plt.figure(figsize=(15,5))
	for key in data.keys():
		tmp = data[key]
		if key == refrun:
			plt.plot(tmp['time'], tmp['uc'], label=tmp['label'], color='k')
		else:
			plt.plot(tmp['time'], tmp['uc'], label=tmp['label'])
	plt.grid()
	plt.title(r'(a) NONES: '+uctitle+r' (m s$^{-1}$)')
	plt.legend()
	plt.savefig(f'Figure_SI_nones_uc')

	# 1. Ice-shelf melting
	plt.figure(figsize=(15,5))
	for key in data.keys():
		tmp = data[key]
		if key == refrun:
			plt.plot(tmp['time'], tmp['isfw'], label=tmp['label'], color='k')
		else:
			plt.plot(tmp['time'], tmp['isfw'], label=tmp['label'])
	plt.grid()
	plt.title(r'(b) NONES: '+isfwtitle+r' (Gt yr$^{-1}$)')
	plt.legend()
	plt.savefig(f'Figure_SI_nones_isfw')

	#==

	quit()

#==

