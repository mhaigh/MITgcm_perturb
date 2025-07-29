# NAT_FIGS.py

import os

import netCDF4 as nc
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

DATA_ROOT = '/data/oceans_output/shelf/michai/mitgcm/'
NC_DATA_ROOT = 'data/'
GRID_PATH = 'data/grid/'

black = 'k'; blue = 'cornflowerblue'; red = 'indianred'; green = 'seagreen'; grey = 'gray'
cmapMean = 'plasma'; cmapAnom = 'coolwarm'; white = 'w'; cyan = 'cyan'; darkgrey = 'dimgrey'
#blue = 'royalblue'; red = 'darkred'; green = 'darkgreen'


# Time constants
secs_per_year = 86400. * 365.25
secs_per_month = secs_per_year / 12.

# Factor for converting kg/m2/s to m/yr.
fwConversion = 1.e-3 * secs_per_year

# Lats for Amundsen Sea subregion.
EASlons = [235, 262]; EASlats = [-75.5, -70]

LATS = [-75.5, -68]
LONS = [225, 262]

# Area of integration. Set longitudinal ranges on/off the shelf.
#lonsN=[225., 245]; lonsS=[233,258.5]; Nrange=58
lonsN=[230., 245]; lonsS=[233,258.5]; Nrange=58
		
#==========================================================

ttt = True
fff = False

#==

# Figure 1
FIGURE1 = fff
if FIGURE1:
	
	# Data needed in this block:
	# Uvel & Vvel, full nc files or npy files at specific level.
	# Theta, full nc file of npy file at specific level.
	# EXFuwind, nc file or npy file.
	# SIfw, nc file or npy file.
	# ISfw, nc file or npy file.
	
	level = 16
	run = 'PAS_668'
	
	# Colorbar params
	vminFW = -10.; vmaxFW = -vminFW
	vminT = -2.; vmaxT = -vminT
	
	pathr = DATA_ROOT + run + '/run/'
	pathp = DATA_ROOT + run + '/post/'
	grid = Grid_PAS(GRID_PATH)	
	bathy = grid.bathy
	draft = grid.draft
	
	iceC = grid.iceC
	iceC = np.where(iceC==0, np.nan, iceC)

	X = grid.XC; Y = grid.YC
	Xl = X[0,:]; Yl = Y[:,0]
	Ny2 = len(Yl)//2; Nx2 = len(Xl)//2
	lat_0=Yl[Ny2]; lon_0=Xl[Nx2]

	lats = LATS; lons = LONS
	latsi = grid.getIndexFromLat(lats); lonsi = grid.getIndexFromLon(lons)

	X, Y = grid.XYsubr(lons, lats)
	Z = grid.RC.squeeze()
	
	##
	# FIG 1A DATA.
	
	# First check if fully processed nc files exist.
	# If one does, assume others do too.
	Tnc = 'time_mean_ocean_potential_temperature_at_455m_depth.nc' 
	if os.path.isfile(f'{NC_DATA_ROOT}Figure_1/{Tnc}'):
		T = nc.Dataset(f'{NC_DATA_ROOT}Figure_1/{Tnc}')
		T = T['theta'][:]
		u = nc.Dataset(f'{NC_DATA_ROOT}Figure_1/time_mean_zonal_ocean_velocity_at_455m_depth.nc')
		u = u['u'][:]
		v = nc.Dataset(f'{NC_DATA_ROOT}Figure_1/time_mean_meridional_ocean_velocity_at_455m_depth.nc')
		v = v['v'][:]
		uw = nc.Dataset(f'{NC_DATA_ROOT}Figure_1/time_mean_2m_atmospheric_zonal_wind.nc')
		uw = uw['uw'][:]
		vw = nc.Dataset(f'{NC_DATA_ROOT}Figure_1/time_mean_2m_atmospheric_meridional_wind.nc')
		vw = vw['vw'][:]
		fw = nc.Dataset(f'{NC_DATA_ROOT}Figure_1/time_mean_sea_ice_and_ice_shelf_freshwater_flux.nc')
		fw = fw['fw'][:]
		
	# Else numpy and un-processed model output require processing.
	else:
		# Load velocities & temp. Take time mean.
		# UVEL
		if os.path.isfile(f'{pathp}UVEL_z{level}.npy'):
			u = np.load(f'{pathp}UVEL_z{level}.npy').mean(axis=0)
		elif os.path.isfile(f'{pathr}stateUvel.nc'):
			u = readVariable('UVEL', pathr, file_format='nc', meta=False, zz=level).mean(axis=0)
		else:
			print('Missing UVEL file. Quitting'); quit()
		# VVEL
		if os.path.isfile(f'{pathp}VVEL_z{level}.npy'):
			v = np.load(f'{pathp}VVEL_z{level}.npy').mean(axis=0)
		elif os.path.isfile(f'{pathr}stateVvel.nc'):
			v = readVariable('VVEL', pathr, file_format='nc', meta=False, zz=level).mean(axis=0)
		else:
			print('Missing VVEL file. Quitting'); quit()
		# THETA
		if os.path.isfile(f'{pathp}THETA_z{level}.npy'):
			T = np.load(f'{pathp}THETA_z{level}.npy').mean(axis=0)
		elif os.path.isfile(f'{pathr}stateTheta.nc'):
			T = readVariable('THETA', pathr, file_format='nc', meta=False, zz=level).mean(axis=0)
		else:
			print('Missing THETA file. Quitting'); quit()
		
		# Process data	

		T = tools.boundData(T, vminT, vmaxT, 0.9999)
		T[-10,0] = vminT; T[-11,1] = vmaxT
		
		u = tools.interp(u, 'u'); v = tools.interp(v, 'v')
		lim = 0.1
		u = tools.boundData(u, -lim, lim, 0.9999); v = tools.boundData(v, -lim, lim, 0.9999)

		##
		# FIG 1B DATA
		# Load velocities & temp. Take time mean.
		# UWIND
		if os.path.isfile(f'{pathp}EXFuwind.npy'):
			uw = np.load(f'{pathp}EXFuwind.npy').mean(axis=0)
		elif os.path.isfile(f'{pathr}stateExf.nc'):
			uw = readVariable('EXFuwind', pathr, file_format='nc', meta=False, var2D=True).mean(axis=0)
		else:
			print('Missing UWIND file. Quitting'); quit()
		# VWIND
		if os.path.isfile(f'{pathp}EXFvwind.npy'):
			vw = np.load(f'{pathp}EXFvwind.npy').mean(axis=0)
		elif os.path.isfile(f'{pathr}stateExf.nc'):
			vw = readVariable('EXFvwind', pathr, file_format='nc', meta=False, var2D=True).mean(axis=0)
		else:
			print('Missing VWIND file. Quitting'); quit()
		# SIFW
		if os.path.isfile(f'{pathp}SIfw.npy'):
			fw = np.load(f'{pathp}SIfw.npy').mean(axis=0)
		elif os.path.isfile(f'{pathr}state2D.nc'):
			fw = readVariable('SIfwmelt', pathr, file_format='nc', meta=False, var2D=True).mean(axis=0)
			fw += readVariable('SIfwfrz', pathr, file_format='nc', meta=False, var2D=True).mean(axis=0)
		else:
			print('Missing SIFW file. Quitting'); quit()
		# ISFW
		if os.path.isfile(f'{pathp}SHIfwFlx.npy'):
			fw -= np.load(f'{pathp}SHIfwFlx.npy').mean(axis=0)
		elif os.path.isfile(f'{pathr}state2D.nc'):
			fw -= readVariable('SHIfwFlx', pathr, file_format='nc', meta=False, var2D=True).mean(axis=0)
		else:
			print('Missing ISFW file. Quitting'); quit()
		fw *= fwConversion
		
		# Process data	
		fw = tools.boundData(fw, vminFW, vmaxFW, 0.9999)

		## END READING DATA
			
		# Get subregion and mask

		u = tools.getSubregionXY(u, latsi, lonsi)
		v = tools.getSubregionXY(v, latsi, lonsi)
		T = tools.getSubregionXY(T, latsi, lonsi)
		uw = tools.getSubregionXY(uw, latsi, lonsi)
		vw = tools.getSubregionXY(vw, latsi, lonsi)
		fw = tools.getSubregionXY(fw, latsi, lonsi)

		# Mask draft
		u = ptt.maskDraftXY(u, grid, 0, timeDep=False, subregion=True, lats=latsi, lons=lonsi)
		v = ptt.maskDraftXY(v, grid, 0, timeDep=False, subregion=True, lats=latsi, lons=lonsi)
		T = ptt.maskDraftXY(T, grid, 0, timeDep=False, subregion=True, lats=latsi, lons=lonsi)
		uw = ptt.maskDraftXY(uw, grid, 0, timeDep=False, subregion=True, lats=latsi, lons=lonsi)
		vw = ptt.maskDraftXY(vw, grid, 0, timeDep=False, subregion=True, lats=latsi, lons=lonsi)
		
		# Mask bathy
		u = ptt.maskBathyXY(u, grid, level, timeDep=False, subregion=True, lats=latsi, lons=lonsi)
		v = ptt.maskBathyXY(v, grid, level, timeDep=False, subregion=True, lats=latsi, lons=lonsi)
		T = ptt.maskBathyXY(T, grid, 0, timeDep=False, subregion=True, lats=latsi, lons=lonsi)
		T = ptt.maskBathyXY(T, grid, level, timeDep=False, subregion=True, lats=latsi, lons=lonsi)
		uw = ptt.maskBathyXY(uw, grid, 0, timeDep=False, subregion=True, lats=latsi, lons=lonsi)
		vw = ptt.maskBathyXY(vw, grid, 0, timeDep=False, subregion=True, lats=latsi, lons=lonsi)
		fw = ptt.maskBathyXY(fw, grid, 0, timeDep=False, subregion=True, lats=latsi, lons=lonsi)

	# End processing data.
	
	#==

	# Grid-related variables
	bathy = tools.getSubregionXY(bathy, latsi, lonsi)
	draft = tools.getSubregionXY(draft, latsi, lonsi)

	iceC = tools.getSubregionXY(iceC, latsi, lonsi)
	land = np.where(bathy<0,0,1)
	
	# This creates land mask with appropriate shade of grey masking.
	nY, nX = bathy.shape
	land = np.zeros((nY, nX, 4))
	for j in range(nY):	
		for i in range(nX):
			if bathy[j,i] == 0:
				land[j,i,:3] = 0.4
			else:
				land[j,i,3] = 1.
				
	iceC = ptt.maskBathyXY(iceC, grid, 0, timeDep=False, subregion=True, lats=latsi, lons=lonsi)
	bathy = ptt.maskBathyXY(bathy, grid, level, timeDep=False, subregion=True, lats=latsi, lons=lonsi)
	bathy = ptt.maskDraftXY(bathy, grid, level, timeDep=False, subregion=True, lats=latsi, lons=lonsi)
	bathyC = ptt.makeBathyContour(bathy, X, Y)

	iceMask = np.ma.getmask(iceC)
	Tmask = np.ma.getmask(bathy)

	
	##
	# PREPARE ALL PLOTTING DATA
	
	# Sample rate
	d = 8
	Xd = X[::d, ::d]; Yd = Y[::d, ::d]
	
	paras = [[-75,  -73,  -71, -69]]*2
	merids = [[230, 240, 250, 260]]*2
	gridLinesWidth = 0	

	# List for ice shelves
	#r = 1.1e6
	r = 2.2e6
	#xloc1 = [2.66e5, 1.55e6+r, 0.9e6+r, 4.6e5+r, 1.85e5+r, 1.52e6+r, 1.4e6+r]
	xloc1 = [9.66e5, 1.52e6+r, 0.9e6+r, 3.8e5+r, 1.6e5+r, 1.52e6+r, 1.4e6+r]
	yloc1 = [4.e5, 1.7e5, 1.6e5, 1.7e5, 3.85e5, 8.e5, 1.2e6]
	text1 = ['GETZ', 'PIG', 'THW', 'CRO', 'DOT', 'COS', 'ABB']
	fontdictIce = {'fontsize':15, 'color':'w'}
	fontdict1 = [fontdictIce]*7
	
	# List for contour labels
	xloc2 = []
	yloc2 = []
	text2 = []
	fontdict2 = []
	
	text0 = {'text':text1+text2, 'xloc':xloc1+xloc2, 'yloc':yloc1+yloc2, 'fontdict':fontdict1+fontdict2}
	text_data = text0
	
	lw = 2.5
	xw = -125; ysw = -72.85; ynw = -71.85
	xe = -108; yse = -71.6; yne = -70.6
	plotvlines = [(xw, ysw, ynw, 'blue', lw), (xe, yse, yne, 'blue', lw)]
	
	r2 = 1.e6
	dxArrow = 0.8e5; dxArrowDiag = 5.67e4
	GetzArrow = [.65e6+r2, 1.50e6, dxArrowDiag, -dxArrowDiag, 'r']
	arrows = [[1.58e6+r2, 8.2e5, dxArrow, 0, 'r'], [.93e6+r2, 1.28e6, -dxArrowDiag, dxArrowDiag, 'b'], [1.18e6+r2, 1.55e6, dxArrow, 0, 'b']]

	u = [u[::d, ::d], uw[::d, ::d]]
	v = [v[::d, ::d], vw[::d, ::d]]
	Xd = [Xd, Xd]
	Yd = [Yd, Yd]
	contourf = [T, fw]
	contourfNlevels = [11, 11]
	cbarTicks = [[-2,-1,0,1,2], [-10,-5,0,5,10]]
	XX = [X, X]
	YY = [Y, Y]
	
	contour = [bathyC, bathyC]
	contourLevels = [[-1000], [-1000]]
	cc = [cyan, cyan]
	lw = [1.5, 1.5]
	
	contour2 = [bathy, draft]
	contourLevels2=[[-500],[-1]]
	cc2 = [white, black]
	lw2 = [1.5, 3]
	
	# Don't need contour outlining Marie Byrd plateau anymore.
	contour3 = [None, None]
	contourLevels3 = [[], [-3000]]
	cc3 = ['', black]
	lw3 = [1., 1.5]
		
	land = [land, land]
	isf = [iceC, None]
	
	arrows = [arrows, None]
	text_data = [text_data, None]
	plotvlines = [plotvlines, None]
	
	titles = [r'(a) Time-mean flow & pot. temp. ($^{\circ}$C) at Z = -455 m', \
	r'(b) Time-mean winds & FW fluxes (m yr$^{-1}$)']
	fstitle = 18
	fontsize = 16
	
	qcolor = [black, white]
	qunits = ['m/s']*2
	scale = [1., 100.]
	qx=0.3; qy=0.03;
	qs = [0.05, 4] 

	width = [0.003]*2
	headwidth = [5.]*2
	headlength = [5.]*2
	headaxislength = [3]*2
	qcolor = [black, white]
	qlabelx = [0.3]*2
	qlabely = [0.028]*2
	
	vmin = [vminT, vminFW]
	vmax = [vmaxT, vmaxFW]

	#==

	pt.quiver1byNBasemap(u, v, Xd, Yd, lat_0, lon_0, contourf=contourf, contourfNlevels=contourfNlevels, cbarTicks=cbarTicks, X=XX, Y=YY, mesh=False, isf=isf, contour=contour, contourLevels=contourLevels, cc=cc, lw=lw, contour2=contour2, contourLevels2=contourLevels2, cc2=cc2, lw2=lw2, contour3=contour3, contourLevels3=contourLevels3, cc3=cc3, lw3=lw3, land=land, titles=titles, fstitle=fstitle, cmap=cmapMean, vmin=vmin, vmax=vmax, parallels=paras, meridians=merids, gridLinesWidth=gridLinesWidth, text_data=text_data,  grid=False, scale=scale, qs=qs, qunits=qunits, extend=['both']*2, width=width, headwidth=headwidth, headlength=headlength, headaxislength=headaxislength, qcolor=qcolor, qlabelx=qlabelx, qlabely=qlabely, figsize=(9.5,11), save=True, outpath='', outname='Figure_1.png', show=True, dpi=200, AntarcticInsetData=[True,None], arrows=arrows, plotvlines=plotvlines, maskColour='0.8')

	quit()

#==
	
# Figure 2
FIGURE2 = fff
if FIGURE2:
	
	alltitle = '(a) ALL simulation timeseries' 	

	ucfname = 'slope_uv_zmax.npy'
	uctitle = r'(b) Undercurrent speed (cm s$^{-1}$)'
	uclim = 0.8
	
	bctitle = r'(c) Baroclinicity (cm s$^{-1}$)'
	bclim = 1.8

	isfwfname = 'isfwSum.npy'	
	isfwtitle = r'(d) Ice-shelf melt (Gt yr$^{-1}$)'
	isfwlim = 90.
	
	sifwfnames = ['']

	runs = [('PAS_668', 'ALL', 'k'), \
		('PAS_668_91-92_winds', 'WINDS', blue), \
		('PAS_668_91-92_thermo', 'THERMO', red), \
		('PAS_668_91-92_none', 'NONE', green), \
		('PAS_668_noFeedback', 'ALL NO FEEDBACK', grey) \
		]

	# Additional params. Running-mean window, undercurrent, sections, sim start year.	
	nn = 48
	sections = ['westGetz', 'westPITW', 'westPITE']
	year0 = 1979
	
	#==
	
	path = DATA_ROOT + runs[0][0] + '/run/'
	grid = Grid_PAS(DATA_ROOT + runs[0][0] + '/run/')
	bathy = grid.bathy

	# Area arrays if any timeseries require area averaging.
	area = grid.RAC * grid.hFacW[0] * grid.hFacS[0]
	areaS, areaN = ptools.maskSouthNorthSlope(area, grid, lonsN=lonsN, lonsS=lonsS, Nrange=Nrange)
	areaNplot = np.where(np.ma.getmask(areaN), 0, areaN)
	areaSplot = np.where(np.ma.getmask(areaS), 0, areaS)
	
	# For each run, load in desired field (e.g., FW or SIheff), average over predefined area.
	data = {}

	for run in runs:
	
		# First try find post-processed data
		ucnc = f'{run[1]}_timeseries_of_along_continental_slope_undercurrent_speed.nc'
		scnc = f'{run[1]}_timeseries_of_along_continental_slope_surface_flow_speed.nc'
		bcnc = f'{run[1]}_timeseries_of_along_continental_slope_baroclinicity.nc'
		windsnc = f'{run[1]}_timeseries_of_along_continental_slope_wind_speed.nc'
		isfwnc = f'{run[1]}_timeseries_of_area_integrated_ice_shelf_melting.nc'
		if os.path.isfile(f'{NC_DATA_ROOT}Figure_2/{ucnc}'):
			uc = nc.Dataset(f'{NC_DATA_ROOT}Figure_2/{ucnc}')
			year = uc['year'][:]
			uc = 100*uc['uc'][:]
			bc = nc.Dataset(f'{NC_DATA_ROOT}Figure_2/{bcnc}')
			bc = 100*bc['bc'][:]
			isfw = nc.Dataset(f'{NC_DATA_ROOT}Figure_2/{isfwnc}')
			isfw = -1.e-12 * isfw['isfw'][:]
			if run[1] == 'ALL':
				winds = nc.Dataset(f'{NC_DATA_ROOT}Figure_2/{windsnc}')
				winds = winds['wind_speed'][:]
				sc = nc.Dataset(f'{NC_DATA_ROOT}Figure_2/{scnc}')
				sc = 100*sc['sc'][:]
			else:
				sc = None
				winds = None
				
		# Else load from un-processed files.
		else:
			# Start with grid and time data
			path = DATA_ROOT + run[0] + '/post/'
			t = np.load(DATA_ROOT + run[0] + '/post/PAS_time.npy')
			year = np.array(ptools.getDecimalTime(t))
			year = year - year[0] + year0
			year = ptools.windowAv(year, n=nn, av=False)
			
			# Undercurrent in cm/s.
			uc = 100*ptools.getUndercurrent(sections=sections, DETREND=True, AV=True, nn=nn, path=path, fname=ucfname)
			
			# Baroclinicity in cm/s.
			bc = 100*ptools.getBarocl(DETREND=True, AV=True, nn=nn, path=path)
			
			# Ice-shelf melting in Gt/yr (+ve ito ocean).
			isfw = -1.e-12 * np.load(path + isfwfname) 
			isfw = ptools.detrend(isfw)
			isfw = ptools.windowAv(isfw, n=nn)
		
			# For ALL, also get winds and surface current.
			if run[1] == 'ALL':
				sc = 100*ptools.getSurfCurrent(sections=sections, DETREND=True, AV=True, nn=nn, path=path)		
				winds = ptools.getSlopeWinds(sections=sections, DETREND=True, AV=True, nn=nn, path=path)
			else:
				sc = None
				winds = None	
		# End collecting data
		
		# Lump all data in dictionary.
		data[run[0]] = {'uc':uc, 'bc':bc, 'isfw':isfw, 'sc':sc, 'winds':winds, 'time':year, 'label':run[1], 'color':run[2]}
	
	#==
	
	fs = 16
	lw = 2.3
	ncol = 5
	
	left = 0.04
	top = 0.95
	bottom = 0.02
	wspace=0.01
	hspace=0.2

	# What to plot for panel 1.
	panel1 = [('uc', 1, red, 'Undercurrent (cm s$^{-1}$)'),\
		  ('bc', 1, green, 'Baroclinicity (cm s$^{-1}$)'),\
		  ('sc', 1, blue, 'Surf. flow (cm s$^{-1}$)'),\
		  ('isfw', 0.01, black, 'Ice-shelf melt (100 Gt yr$^{-1}$)'),\
		  ('winds', 1, grey, 'Winds (m s$^{-1}$)')]
	 
	# PLOT
	plt.figure(figsize=(15,16))
	
	plt.subplot(411)
	tmp = data['PAS_668']
	for var in panel1:	
		plt.plot(tmp['time'], var[1]*tmp[var[0]], label=var[3], color=var[2], linewidth=lw)
	plt.xlim(min(tmp['time']), max(tmp['time']))
	plt.ylim(-1.9, 1.9)
	plt.xticks(fontsize=fs-2)
	plt.yticks(fontsize=fs-2)
	plt.grid()
	plt.title(alltitle, fontsize=fs+2)
	plt.legend(prop={'size':fs-3}, ncol=ncol)

	# Panel 1: undercurrent
	plt.subplot(412)
	for key in data.keys():
		tmp = data[key]
		plt.plot(tmp['time'], tmp['uc'], label=tmp['label'], color=tmp['color'], linewidth=lw)
	plt.xlim(min(tmp['time']), max(tmp['time']))
	plt.ylim(-uclim,uclim)
	plt.xticks(fontsize=fs-2)
	plt.yticks(fontsize=fs-2)
	plt.grid()
	plt.title(uctitle, fontsize=fs+2)
	plt.legend(prop={'size':fs-2}, ncol=ncol)
	
        # Panel 2: baroclinicity
	plt.subplot(413)
	for key in data.keys():
		tmp = data[key]
		if tmp['label'] == 'ALL':
			std = np.mean((tmp['bc']-np.mean(tmp['bc']))**2)**0.5
			plt.axhline(std, color='k', linestyle='--', linewidth=1.2)
			plt.axhline(-std, color='k', linestyle='--', linewidth=1.2)
		plt.plot(tmp['time'], tmp['bc'], label=tmp['label'], color=tmp['color'], linewidth=lw)
	plt.xlim(min(tmp['time']), max(tmp['time']))
	plt.ylim(-bclim, bclim)
	plt.xticks(fontsize=fs-2)
	plt.yticks(fontsize=fs-2)
	plt.grid()
	plt.title(bctitle,fontsize=fs+2)
	plt.legend(prop={'size':fs-2}, ncol=ncol)

	# Panel 3: ice-shelf melting
	plt.subplot(414)
	for key in data.keys():
		tmp = data[key]
		plt.plot(tmp['time'], tmp['isfw'], label=tmp['label'], color=tmp['color'], linewidth=lw)
	plt.xlim(min(tmp['time']), max(tmp['time']))
	plt.ylim(-isfwlim, isfwlim)	
	plt.xticks(fontsize=fs-2)
	plt.yticks(fontsize=fs-2)
	plt.grid()
	plt.title(isfwtitle, fontsize=fs+2)
	plt.legend(prop={'size':fs-2}, ncol=ncol)
	
	# Panel 3: sea-ice freshwater flux south
	#plt.subplot(414)
	#for key in data.keys():
	#	tmp = data[key]
	#	plt.plot(tmp['time'], tmp['sifwS'], label=tmp['label'], color=tmp['color'])
	#plt.ylim(-sifwlim, sifwlim)	
	#plt.grid()
	#plt.title(sifwStitle)
	#plt.legend()
	
	#plt.tight_layout()
	plt.subplots_adjust(left=left, right=1-left, bottom=bottom, top=top, wspace=wspace, hspace=hspace)	

	plt.savefig('Figure_2')
	plt.show()

	quit()
	
#==

FIGURE3 = fff
if FIGURE3:

	##
	# OPTIONS
	VAR = ['slopeRho', 1, 'slope Rho', (-0.05, 0.05)]
	##

	runs = [('PAS_668', 'ALL', [(0,0), (0,1)]), \
		('PAS_668_91-92_winds', 'WINDS', [(0,2), (0,3)]), \
		('PAS_668_91-92_thermo', 'THERMO', [(1,2), (1,3)]), \
		('PAS_668_91-92_none', 'NONE', [(2,2), (2,3)]), \
		('PAS_668_91-92_noFeedback', 'ALL NO FEEDBACK', [(1,0), (1,1)])]
	
	#==
	
	allSections = ['westGetz','Getz','westPITW','PITW','westPITE']
	allSection_indices = ptools.getSectionIndices(allSections)

	# Start and end times of desired ranges.
	opt = 'range'; thresh = 1
	s1 = [20, 40]; e1 = [60, 100]
	s2 = [100, 160]; e2 = [160, 200]
	s3 = [280, 320]; e3 = [370, 390]
	s4 = [390, 410]; e4 = [420, 460]

	letters = [['(a)', '(b)', '(c)', '(d)'], \
	           ['(e)', '(f)', '(g)', '(h)'], \
		   ['(i)', '(j)', '(k)', '(l)']]

	SUM_COMPS = True
	DETREND = True
	nn = 48
	
	year0 = 1979

	cp = 3900; rho0 = 1028.5
	lhf = 334.e3
	dt_month = 86400. * 365.25 / 12
	dt_year = 86400. * 365.25
	offsetStr = ''
	vmin = VAR[3][0]; vmax = VAR[3][1]

	#==
	
	pathpref =  DATA_ROOT + 'PAS_668/post/'	
	pathrref = DATA_ROOT + runs[0][0] + '/run/'
	grid = Grid_PAS(pathrref)
	slopeBathy = np.load(pathpref+'slopeBathy.npy')
	
	yslim = 50; ynlim = 50; nY = ynlim + yslim
	dY = 3.6e0; Y = np.linspace(-yslim*dY, ynlim*dY, nY)
	zs = [0, -1000]; zz = grid.getIndexFromDepth(zs)
	Z = grid.RC.squeeze()[zz[0]:zz[1]]
	nZ = Z.shape[0]	

	# Load time
	t = np.load(pathpref + 'PAS_time.npy')
	year = np.array(ptools.getDecimalTime(t))
	year = year - year[0] + year0
	year = ptools.windowAv(year, n=nn, av=False)

	# Get indices for starts/ends of averaging periods
	#uc = ptools.getUndercurrent(DETREND=True, AV=True, nn=nn, path=pathpref, fname='slope_uv_zmax.npy')	
	bc = ptools.getBarocl(DETREND=True, AV=True, nn=nn, path=pathpref)

	# Get ranges for computing composites. If unsure check outputted timeseries plot.
	# This one for composites taken between +/- thresh std
	s1i, e1i, s2i, e2i, s3i, e3i, s4i, e4i = ptools.computeFourCompositeRanges(s1, e1, s2, e2, s3, e3, s4, e4, thresh=thresh, year=year, tseries=bc, opt=opt, plotCheck=False)

	#==

	# Initialise list this way so we can overwrite.
	COMPS = [[np.zeros((nZ,nY)) for i in range(4)] for j in range(3)]
	UCOMPS = [[np.zeros((nZ,nY)) for i in range(4)] for j in range(3)]
	insetTitle = [[None for i in range(4)] for j in range(3)]
	feedbackAnomPos = np.zeros((nZ,nY))
	feedbackAnomNeg = np.zeros((nZ,nY))
	UfeedbackAnomPos = np.zeros((nZ,nY))
	UfeedbackAnomNeg = np.zeros((nZ,nY))

	# Start main loop through each run
	for run in runs:
	
		print(run[0])
		pathr = DATA_ROOT + run[0] + '/run/'
		pathp = DATA_ROOT + run[0] + '/post/'
	
		# If processed nc files exist...
		rhonc = f'{run[1]}_cross_slope_density_anomaly_composites.nc'
		unc = f'{run[1]}_along_slope_velocity_anomaly_composites.nc'
		if os.path.isfile(f'{NC_DATA_ROOT}Figure_3/{unc}'):
			print(f'Reading {run[1]} netCDF files')
			u = nc.Dataset(f'{NC_DATA_ROOT}Figure_3/{unc}')
			ucomp1 = u['u_slope_high_baroclinicity'][:]
			ucomp2 = u['u_slope_low_baroclinicity'][:]
			rho = nc.Dataset(f'{NC_DATA_ROOT}Figure_3/{rhonc}')
			comp1 = rho['rho_high_baroclinicity'][:]
			comp2 = rho['rho_low_baroclinicity'][:]
					
		# Else load un-processed data
		else:
			# LOAD DATA
			data = np.load(pathp+'slopeRho.npy')[...,allSection_indices].mean(axis=-1)
			udata = np.load(pathp+'slope_uv_xyz.npy')[...,allSection_indices].mean(axis=-1)
		
			#==
			
			# TREAT DATA & GET COMPOSITES
			data = ptools.demean(data)
			udata = ptools.demean(udata)
			if DETREND:
				data = ptools.detrend(data)
				udata = ptools.detrend(udata)#, interceptFlag=0)
			data = ptools.windowAv(data, n=nn)
			udata = ptools.windowAv(udata, n=nn)

			# Get composite
			comp1 = np.mean(data[s1i:e1i], axis=0)
			comp2 = np.mean(data[s2i:e2i], axis=0)
			comp3 = np.mean(data[s3i:e3i], axis=0)
			comp4 = np.mean(data[s4i:e4i], axis=0)

			ucomp1 = np.mean(udata[s1i:e1i], axis=0)	
			ucomp2 = np.mean(udata[s2i:e2i], axis=0)
			ucomp3 = np.mean(udata[s3i:e3i], axis=0)
			ucomp4 = np.mean(udata[s4i:e4i], axis=0)

			    # This replaces comp1 and comp2 with average over both instances
		            # of fast/slow undercurrent. Comment out as default 
			if SUM_COMPS:
				comp1 = np.sum(data[s1i:e1i],axis=0)+np.sum(data[s3i:e3i], axis=0)
				comp1 = comp1 / (e1i-s1i + e3i-s3i)
				comp2 = np.sum(data[s2i:e2i],axis=0)+np.sum(data[s4i:e4i], axis=0)
				comp2 = comp2 / (e2i-s2i + e4i-s4i)
				ucomp1 = np.sum(udata[s1i:e1i],axis=0)+np.sum(udata[s3i:e3i], axis=0)
				ucomp1 = ucomp1 / (e1i-s1i + e3i-s3i)
				ucomp2 = np.sum(udata[s2i:e2i],axis=0)+np.sum(udata[s4i:e4i], axis=0)
				ucomp2 = ucomp2 / (e2i-s2i + e4i-s4i)

			# Set bounds
			comp1 = tools.boundData(comp1, vmin, vmax, d=1.e-9)
			comp2 = tools.boundData(comp2, vmin, vmax, d=1.e-9)
			comp3 = tools.boundData(comp3, vmin, vmax, d=1.e-9)
			comp4 = tools.boundData(comp4, vmin, vmax, d=1.e-9)
			
			# Mask
			bathyMask = np.zeros(comp1.shape)
			for yi in range(nY):
				bathyMask[:, yi] = np.max(slopeBathy[yi, allSection_indices], axis=-1)
			for zi in range(len(Z)):
				bathyMask[zi,:] -= Z[zi]	

			comp1 = np.ma.masked_where(bathyMask>0, comp1)
			comp2 = np.ma.masked_where(bathyMask>0, comp2)
			comp3 = np.ma.masked_where(bathyMask>0, comp3)
			comp4 = np.ma.masked_where(bathyMask>0, comp4)

			ucomp1 = np.ma.masked_where(bathyMask>0, ucomp1)
			ucomp2 = np.ma.masked_where(bathyMask>0, ucomp2)
			ucomp3 = np.ma.masked_where(bathyMask>0, ucomp3)
			ucomp4 = np.ma.masked_where(bathyMask>0, ucomp4)
		# End processing data

		# Put data into list of lists
		if SUM_COMPS:
			row1 = run[2][0][0]; col1 = run[2][0][1]
			row2 = run[2][1][0]; col2 = run[2][1][1]
			COMPS[row1][col1] = comp1
			COMPS[row2][col2] = comp2
			UCOMPS[row1][col1] = ucomp1
			UCOMPS[row2][col2] = ucomp2
			insetTitle[row1][col1] = f'{letters[row1][col1]} {run[1]}\nmax baroclinicity'
			insetTitle[row2][col2] = f'{letters[row2][col2]} {run[1]}\nmin baroclinicity' 
			if run[1] == 'ALL':
				feedbackAnomPos += comp1			
				feedbackAnomNeg += comp2
				UfeedbackAnomPos += ucomp1
				UfeedbackAnomNeg += ucomp2
			elif run[1] == 'ALL NO FEEDBACK':
				feedbackAnomPos -= comp1
				feedbackAnomNeg -= comp2
				UfeedbackAnomPos -= ucomp1
				UfeedbackAnomNeg -= ucomp2

		# THIS OPTION WON'T WORK WITH PLOTTING OPTS BELOW
		else:
			print('SUM_COMPS=False may lead to errors.')
			COMPS.append([comp1, comp2, comp3, comp4])
	

		# Last, print min/max values
		PRINT_MIN_MAX = True
		if PRINT_MIN_MAX:
			print(f'{run[1]} max / min:')
			print(f'{np.min(comp1)} / {np.max(comp2)}')

	# End  run loop
	#==
	
	# Put NO FEEDBACK anomaly into COMPS
	feedbackAnomPos = np.ma.masked_where(bathyMask>0, feedbackAnomPos) 
	feedbackAnomNeg = np.ma.masked_where(bathyMask>0, feedbackAnomNeg)
	COMPS[2][0] = feedbackAnomPos
	COMPS[2][1] = feedbackAnomNeg
	UfeedbackAnomPos = np.ma.masked_where(bathyMask>0, UfeedbackAnomPos)
	UfeedbackAnomNeg = np.ma.masked_where(bathyMask>0, UfeedbackAnomNeg)
	UCOMPS[2][0] = UfeedbackAnomPos
	UCOMPS[2][1] = UfeedbackAnomNeg
	insetTitle[2][0] = f'{letters[2][0]} Diff. (a)-(e)\nmax baroclinicity'
	insetTitle[2][1] = f'{letters[2][1]} Diff (b)-(f)\nmin baroclinicity'

	#==

	#PLOT

	DASH_NEG_CONTOURS = True
	contourLevels = [-.004,-.003,-.002,-.001,0,.001,.002,.003,.004]
	contourfNlevels = 11

	xlabel = r'$y_{\text{slope}}$ (km)'
	ylabel = 'Depth (m)'
	xticks = None
	xticklabels = None 
	yticks = None
	yticklabels = None

	spbottom = 0.13
	spleft = 0.1
	spright = 1-spleft

	insetTitleX = 0.05
	insetTitleY = 0.17
	suptitle = r'Cross-slope density composites (kg m$^{-3}$)'
	ylim = (-600,0)	
	figsize = (9.5,9.5)

	pt.plotMbyN(COMPS, X=Y, Y=Z, xlabel=xlabel, ylabel=ylabel, ylim=ylim, \
	vmin=vmin, vmax=vmax, contourfNlevels=contourfNlevels, \
	xticks=xticks, yticks=yticks, xticklabels=xticklabels, yticklabels=yticklabels, \
	contourLevels=contourLevels, contour=UCOMPS, DASH_NEG_CONTOURS=DASH_NEG_CONTOURS, \
	suptitle=suptitle, insetTitle=insetTitle, \
	insetTitleX=insetTitleX, insetTitleY=insetTitleY, \
	spbottom=spbottom, spleft=spleft, spright=spright, figsize=figsize, \
	show=False, save=True, outname=f'Figure_3.png')

	quit()

#==

FIGURE4 = ttt
if FIGURE4:

	##
	# OPTIONS
	VAR = [['SIfwfrz','SIfwmelt'],[fwConversion,fwConversion],'SI FW flux', (-.8,.8), 185]

	# Leave as None if no vector wanted.
	# VECVAR = [(UVEC,VVEC),LEVEL(None if 2D),TITLE,(VMIN,VMAX),d,SCALE,QKEY]
	VECVAR = None
	VECVAR = [('EXFuwind', 'EXFvwind'),None,'Winds',(-.2,.2),12,6.,.4]	
	##

	DO_TIMESERIES = True

	runs = [('PAS_668', 'ALL', 'k'), \
		('PAS_668_91-92_winds', 'WINDS', blue), \
		('PAS_668_91-92_thermo', 'THERMO', red), \
		('PAS_668_91-92_none', 'NONE', green), \
		]

	# Subregion
	SUBRGN = True; lonsi=None; latsi=None

	##
	# These for 4-year running mean barocl, thresh=1.
	opt = 'range'; thresh = 1
	s1 = [20, 40]; e1 = [60, 100]
	s2 = [100, 160]; e2 = [160, 200]
	s3 = [280, 320]; e3 = [370, 390]
	s4 = [390, 410]; e4 = [420, 460]
	##

	DETREND = True
	FILTER = False
	REMOVE_AREA_MEAN = False
	OFFSET = True
	# If SUM_COMPS=True, computes composites averaged over both instances of max/min
	SUM_COMPS = True 
	
	nn = 48
	year0 = 1979
	dt_month = 86400. * 365.25 / 12
	dt_year = 86400. * 365.25
	offsetStr = ''

	#==
	
	pathpref =  DATA_ROOT + 'PAS_668/post/'	
	pathrref =  DATA_ROOT + 'PAS_668/run/'
	grid = Grid_PAS(pathrref)
	bathy = grid.bathy
	draft = grid.draft
	X = grid.XC; Y = grid.YC

	# Area arrays for averaging around the coast or shelf.
	lonsArea=[225,260]; Nrange=58
	lonsS=[233,261.5]
	area = grid.RAC * grid.hFacW[0] * grid.hFacS[0]
	areaS, areaN = ptools.maskSouthNorthSlope(area, grid, lonsN=lonsN, lonsS=lonsS, Nrange=Nrange)
	areaShelf = np.where(np.ma.getmask(areaS), 0, areaS)
	coast = ptools.coastalMask(grid, maskRad=15, omitBurkeIsland=True)
	areaCoast = coast*area
	areaCoast = np.where(X<lonsArea[0], 0, areaCoast)
	areaCoast = np.where(X>lonsArea[1], 0, areaCoast)

	mask = np.ma.getmask(areaS)
	areaS = np.where(mask,1,0)

	# Ice shelf draft
	iceC = grid.iceC
	iceC = np.where(iceC==0, np.nan, iceC)

	# Update LATS
	LATS = [-75.5, -70]
	if SUBRGN:
		latsi = grid.getIndexFromLat(LATS); lonsi = grid.getIndexFromLon(LONS)
		bathy = tools.getSubregionXY(bathy, latsi, lonsi)
		draft = tools.getSubregionXY(draft, latsi, lonsi)
		iceC = tools.getSubregionXY(iceC, latsi, lonsi)
		X, Y = grid.XYsubr(LONS, LATS)

	# Data vmin/vmax
	vmin = VAR[3][0]; vmax = VAR[3][1]

	#==			

	# Prep lists for run loop
	
	COMPS = []
	if VECVAR is None:
		UCOMPS = [[None for i in range(2)] for j in range(4)]
		VCOMPS = [[None for i in range(2)] for j in range(4)]
		Xd = None; Yd = None
		qs = None; scale = None
	else:
		uvmin = VECVAR[3][0]; uvmax = VECVAR[3][1]
		d = VECVAR[4]
		scale = VECVAR[5]
		qs = VECVAR[6]
		UCOMPS = []
		VCOMPS = []
		Xd = X[::d, ::d]; Yd = Y[::d, ::d]

	# START RUN LOOP

	for run in runs:
	
		print(run[0])
		pathp = DATA_ROOT + run[0] + '/post/'

		# First check if processed nc file exists
		fwnc = f'{run[1]}_sea_ice_freshwater_flux_anomaly_composites.nc'
		if os.path.isfile(f'{NC_DATA_ROOT}Figure_4/{fwnc}'):	
			print(f'Reading {run[1]} netCDF files')
			fw = nc.Dataset(f'{NC_DATA_ROOT}Figure_4/{fwnc}')
			comp1 = fw['sifw_pre_high_bc'][:]
			comp2 = fw['sifw_pre_low_bc'][:]
			if run[1] == 'ALL':
				uwnc = 'ALL_2m_atmospheric_zonal_wind_anomaly_composites.nc'
				uw = nc.Dataset(f'{NC_DATA_ROOT}Figure_4/{uwnc}')
				ucomp1 = uw['uw_pre_high_bc']
				ucomp2 = uw['uw_pre_low_bc']
				vwnc = 'ALL_2m_atmospheric_meridional_wind_anomaly_composites.nc'
				vw = nc.Dataset(f'{NC_DATA_ROOT}Figure_4/{vwnc}')
				vcomp1 = vw['vw_pre_high_bc']
				vcomp2 = vw['vw_pre_low_bc']
			bc = None # If data loaded from netcdf, baroclinicty not read in yet.
			# End read nc files
			
		# Else read and process data
		else:
		
			# LOAD DATA
			fname = VAR[0]
			if isinstance(fname, str):
				data = VAR[1]*readVar(fname, pathp)
			else:
				scaleFac = ptt.makeList(VAR[1], len(fname))
				data = scaleFac[0]*readVar(fname[0], pathp)
				for fi in range(1,len(fname)):
					data += scaleFac[fi]*readVar(fname[fi], pathp)

			# Load time
			t = np.load(pathp + 'PAS_time.npy')
			year = np. array(ptools.getDecimalTime(t))
			year = year - year[0] + year0
			year = ptools.windowAv(year, n=nn, av=False)

			### This section could be taken out of loop.
			# Get indices for starts/ends of averaging periods
			bc = ptools.getBarocl(DETREND=True, AV=True, nn=nn, path=pathpref)

			# Get ranges for computing composites. If unsure check outputted timeseries plot.
			# This one for composites taken between +/- thresh std
			s1i, e1i, s2i, e2i, s3i, e3i, s4i, e4i = ptools.computeFourCompositeRanges(s1, e1, s2, e2, s3, e3, s4, e4, thresh=thresh, year=year, tseries=bc, opt=opt, plotCheck=True)

			# Look at periods of time leading up to maxima/minima, instead of the maxima/minima.
			if OFFSET:
				offsetStr = '_offset'
				e4i = s4i
				s4i = e3i
				e3i = s3i
				s3i = e2i
				e2i = s2i
				s2i = e1i
				e1i = s1i
				s1i = 0
			###
					
			if VECVAR is not None:
				if VECVAR[1] == None:
					var2D = True
				else:
					var2D = False
				udata = readVar(VECVAR[0][0], pathp, zz=VECVAR[1], var2D=var2D)
				vdata = readVar(VECVAR[0][1], pathp, zz=VECVAR[1], var2D=var2D)
				
			#==
			
			# TREAT DATA & GET COMPOSITES
			data = ptools.demean(data)
			if DETREND:
				data = ptools.detrend(data)#, interceptFlag=0)
			data = ptools.windowAv(data, n=nn)
			
			if VECVAR is not None:
				udata = ptools.demean(udata)
				vdata = ptools.demean(vdata)
				#if DETREND:
				#	udata = ptools.detrend(udata)
				#	vdata = ptools.detrend(vdata)
				udata = ptools.windowAv(udata, n=nn)
				vdata = ptools.windowAv(vdata, n=nn)
			
			if SUBRGN:
				data = tools.getSubregionXY(data, latsi, lonsi)
				if VECVAR is not None:
					udata = tools.getSubregionXY(udata, latsi, lonsi)
					vdata = tools.getSubregionXY(vdata, latsi, lonsi)
					
			# Get composite
			comp1 = np.mean(data[s1i:e1i], axis=0)
			comp2 = np.mean(data[s2i:e2i], axis=0)
			comp3 = np.mean(data[s3i:e3i], axis=0)
			comp4 = np.mean(data[s4i:e4i], axis=0)
			
			# This replaces comp1 and comp2 with average over both instances
			# of fast/slow undercurrent. Comment out as default 
			if SUM_COMPS:
				comp1 = np.sum(data[s1i:e1i],axis=0)+np.sum(data[s3i:e3i], axis=0)
				comp1 = comp1 / (e1i-s1i + e3i-s3i)
				comp2 = np.sum(data[s2i:e2i],axis=0)+np.sum(data[s4i:e4i], axis=0)
				comp2 = comp2 / (e2i-s2i + e4i-s4i)

			if VECVAR is not None:
				ucomp1 = np.mean(udata[s1i:e1i], axis=0)
				ucomp2 = np.mean(udata[s2i:e2i], axis=0)
				ucomp3 = np.mean(udata[s3i:e3i], axis=0)
				ucomp4 = np.mean(udata[s4i:e4i], axis=0)
				vcomp1 = np.mean(vdata[s1i:e1i], axis=0)
				vcomp2 = np.mean(vdata[s2i:e2i], axis=0)
				vcomp3 = np.mean(vdata[s3i:e3i], axis=0)
				vcomp4 = np.mean(vdata[s4i:e4i], axis=0)

				if SUM_COMPS:
					ucomp1 = np.sum(udata[s1i:e1i],axis=0)+np.sum(udata[s3i:e3i], axis=0)
					ucomp1 = ucomp1 / (e1i-s1i + e3i-s3i)
					ucomp2 = np.sum(udata[s2i:e2i],axis=0)+np.sum(udata[s4i:e4i], axis=0)
					ucomp2 = ucomp2 / (e2i-s2i + e4i-s4i)
					vcomp1 = np.sum(vdata[s1i:e1i],axis=0)+np.sum(vdata[s3i:e3i], axis=0)
					vcomp1 = vcomp1 / (e1i-s1i + e3i-s3i)
					vcomp2 = np.sum(vdata[s2i:e2i],axis=0)+np.sum(vdata[s4i:e4i], axis=0)
					vcomp2 = vcomp2 / (e2i-s2i + e4i-s4i)

			if REMOVE_AREA_MEAN:
				comp1 = comp1 - np.mean(comp1)
				comp2 = comp2 - np.mean(comp2)
				comp3 = comp3 - np.mean(comp3)
				comp4 = comp4 - np.mean(comp4)

			if FILTER:
				print('Spatially smoothing data. TURN OFF if not wanted.')
				comp1 = tools.spatialFilter2D(comp1, 15)
				comp2 = tools.spatialFilter2D(comp2, 15)
				comp3 = tools.spatialFilter2D(comp3, 15)
				comp4 = tools.spatialFilter2D(comp4, 15)

			# Set bounds
			comp1 = tools.boundData(comp1, vmin, vmax, d=1.e-9)
			comp2 = tools.boundData(comp2, vmin, vmax, d=1.e-9)
			comp3 = tools.boundData(comp3, vmin, vmax, d=1.e-9)
			comp4 = tools.boundData(comp4, vmin, vmax, d=1.e-9)
			
			if VECVAR is not None:
				ucomp1 = tools.boundData(ucomp1, uvmin, uvmax, d=1.e-9)
				ucomp2 = tools.boundData(ucomp2, uvmin, uvmax, d=1.e-9)
				ucomp3 = tools.boundData(ucomp3, uvmin, uvmax, d=1.e-9)
				ucomp4 = tools.boundData(ucomp4, uvmin, uvmax, d=1.e-9)
				vcomp1 = tools.boundData(vcomp1, uvmin, uvmax, d=1.e-9)
				vcomp2 = tools.boundData(vcomp2, uvmin, uvmax, d=1.e-9)
				vcomp3 = tools.boundData(vcomp3, uvmin, uvmax, d=1.e-9)
				vcomp4 = tools.boundData(vcomp4, uvmin, uvmax, d=1.e-9)
			
			comp1 = ptt.maskBathyXY(comp1,grid,0,timeDep=False,subregion=SUBRGN,lons=lonsi,lats=latsi)
			comp1 = ptt.maskDraftXY(comp1,grid,0,timeDep=False,subregion=SUBRGN,lons=lonsi,lats=latsi)
			comp2 = ptt.maskBathyXY(comp2,grid,0,timeDep=False,subregion=SUBRGN,lons=lonsi,lats=latsi)
			comp2 = ptt.maskDraftXY(comp2,grid,0,timeDep=False,subregion=SUBRGN,lons=lonsi,lats=latsi)
			comp3 = ptt.maskBathyXY(comp3,grid,0,timeDep=False,subregion=SUBRGN,lons=lonsi,lats=latsi)
			comp3 = ptt.maskDraftXY(comp3,grid,0,timeDep=False,subregion=SUBRGN,lons=lonsi,lats=latsi)
			comp4 = ptt.maskBathyXY(comp4,grid,0,timeDep=False,subregion=SUBRGN,lons=lonsi,lats=latsi)
			comp4 = ptt.maskDraftXY(comp4,grid,0,timeDep=False,subregion=SUBRGN,lons=lonsi,lats=latsi)
			
			if VECVAR is not None:
				ucomp1 = ptt.maskBathyXY(ucomp1,grid,0,timeDep=False,subregion=SUBRGN,lons=lonsi,lats=latsi)
				ucomp1 = ptt.maskDraftXY(ucomp1,grid,0,timeDep=False,subregion=SUBRGN,lons=lonsi,lats=latsi)
				ucomp2 = ptt.maskBathyXY(ucomp2,grid,0,timeDep=False,subregion=SUBRGN,lons=lonsi,lats=latsi)
				ucomp2 = ptt.maskDraftXY(ucomp2,grid,0,timeDep=False,subregion=SUBRGN,lons=lonsi,lats=latsi)
				ucomp3 = ptt.maskBathyXY(ucomp3,grid,0,timeDep=False,subregion=SUBRGN,lons=lonsi,lats=latsi)
				ucomp3 = ptt.maskDraftXY(ucomp3,grid,0,timeDep=False,subregion=SUBRGN,lons=lonsi,lats=latsi)
				ucomp4 = ptt.maskBathyXY(ucomp4,grid,0,timeDep=False,subregion=SUBRGN,lons=lonsi,lats=latsi)
				ucomp4 = ptt.maskDraftXY(ucomp4,grid,0,timeDep=False,subregion=SUBRGN,lons=lonsi,lats=latsi)
				vcomp1 = ptt.maskBathyXY(vcomp1,grid,0,timeDep=False,subregion=SUBRGN,lons=lonsi,lats=latsi)
				vcomp1 = ptt.maskDraftXY(vcomp1,grid,0,timeDep=False,subregion=SUBRGN,lons=lonsi,lats=latsi)
				vcomp2 = ptt.maskBathyXY(vcomp2,grid,0,timeDep=False,subregion=SUBRGN,lons=lonsi,lats=latsi)
				vcomp2 = ptt.maskDraftXY(vcomp2,grid,0,timeDep=False,subregion=SUBRGN,lons=lonsi,lats=latsi)
				vcomp3 = ptt.maskBathyXY(vcomp3,grid,0,timeDep=False,subregion=SUBRGN,lons=lonsi,lats=latsi)
				vcomp3 = ptt.maskDraftXY(vcomp3,grid,0,timeDep=False,subregion=SUBRGN,lons=lonsi,lats=latsi)
				vcomp4 = ptt.maskBathyXY(vcomp4,grid,0,timeDep=False,subregion=SUBRGN,lons=lonsi,lats=latsi)
				vcomp4 = ptt.maskDraftXY(vcomp4,grid,0,timeDep=False,subregion=SUBRGN,lons=lonsi,lats=latsi)
	
		# Put data into list.
		if SUM_COMPS:
			COMPS.append([comp1, comp2])
			if VECVAR is not None:
				UCOMPS.append([ucomp1[::d,::d], ucomp2[::d,::d]])
				VCOMPS.append([vcomp1[::d,::d], vcomp2[::d,::d]])
		else:
			COMPS.append([comp1, comp2, comp3, comp4])
			if VECVAR is not None:
				UCOMPS.append([ucomp1[::d,::d], ucomp2[::d,::d], ucomp3[::d,::d], ucomp4[::d,::d]])
				VCOMPS.append([vcomp1[::d,::d], vcomp2[::d,::d], vcomp3[::d,::d], vcomp4[::d,::d]])
		
		# This makes it so vectors only done for first run in runs.	
		VECVAR = None
		UCOMPS.append([None]*4)
		VCOMPS.append([None]*4)

	# END RUN LOOP
	
	#==
	
	ts_data = None
	ts_time = []; ts_lines = []; ts_labels = []
	if DO_TIMESERIES:
		ts_data = []
		print('Doing timeseries...')
		for run in runs:
			print(run[0])
			
			# First check if processed nc file exists.
			fwnc = f'{run[1]}_cumulative_coastal_sea_ice_freshwater_flux.nc'
			if os.path.isfile(f'{NC_DATA_ROOT}Figure_4/{fwnc}'):
				print(f'Reading {run[1]} sea ice timeseries') 
				fw = nc.Dataset(f'{NC_DATA_ROOT}Figure_4/{fwnc}')
				year = fw['year'][:]
				data = fw['cumulative_coastal_sea_ice_fw_flux'][:]
			
			# Else process from model output.
			else:
				pathp = DATA_ROOT + run[0] + '/post/'

				fname = VAR[0]
				if isinstance(fname, str):
					data = readVar(fname, pathp)
				else:
					data = readVar(fname[0], pathp)
					for fi in range(1,len(fname)):
						data += readVar(fname[fi], pathp)

				data = 86400.*365.25*1.e-12/12 * np.sum(data * areaCoast, axis=(1,2))
				data = ptools.detrend(data)
				data = ptools.windowAv(data, n=nn)
			# End process data
	
			ts_data.append(data)	
			ts_time.append(year)
			ts_lines.append((run[2], 'solid'))
			ts_labels.append(run[1])		
		# END SIFW RUN LOOP

		# Also get baroclinicity for comparison
		if bc is None:
			bcnc = f'ALL_timeseries_of_along_continental_slope_baroclinicity.nc'
			bc = nc.Dataset(f'{NC_DATA_ROOT}Figure_2/{bcnc}')
			bc = bc['bc'][:]
			
		ts_data.append(0.9*VAR[4]*bc/np.max(np.abs(bc)))  
		ts_time.append(year)
		ts_lines.append(('k', 'dashed'))
		ts_labels.append('ALL baroclinicity')	
	# END IF DO_TIMESERIES

	#==

	#PLOT
	
	contourfNlevels = 9

	bathyC = ptt.makeBathyContour(bathy, X, Y)
	contour = bathyC
	contourLevels = [-1000]
	
	coastContour = np.where(areaCoast!=0, 1, 0)
	if SUBRGN:
		coastContour = tools.getSubregionXY(coastContour, latsi, lonsi)
		shelfContour = tools.getSubregionXY(areaS, latsi, lonsi)
	shelfContour[0,:] = 0; shelfContour[1,:] = 1; shelfContour[2,:]

	contour2 = [[None for i in range(2)] for j in range(4)]
	contour2Levels = [[None for i in range(2)] for j in range(4)]
	contour2[3][0] = coastContour
	contour2Levels[3][0] = [0,1]

	xticks = [230, 240, 250, 260]
	yticks = [-75, -73, -71, -69]
	xticklabels = [r'130$^{\circ}$W', r'120$^{\circ}$W', r'110$^{\circ}$W', r'100$^{\circ}$W']
	yticklabels = [r'75$^{\circ}$S', r'73$^{\circ}$S', r'71$^{\circ}$S', r'69$^{\circ}$S']
	extend = 'neither'

	suptitle = r'Sea ice freshwater flux composites (m yr$^{-1}$) & integrated timeseries (Gt)'
	ts_insetTitle = '(i) Cumulative coastal sea-ice freshwater flux (Gt)'

	insetTitle = [[f'(a) {runs[0][1]} pre max baroclinicity', f' (b) {runs[0][1]} pre min baroclinicity'], [f'(c) {runs[1][1]} pre max baroclinicity', f' (d) {runs[1][1]} pre min baroclinicity'], [f'(e) {runs[2][1]} pre max baroclinicity', f' (f) {runs[2][1]} pre min baroclinicity'], [f'(g) {runs[3][1]} pre max baroclinicity', f' (h) {runs[3][1]} pre min baroclinicity']]
	insetTitleX = 0.02

	# General timeseries params. Might not be used.
	ts_xlim = (year[0], year[-1])
	ts_ylim = (-VAR[4], VAR[4])
	if DO_TIMESERIES: 
		hspace=0.05; wspace=0.03; spleft=0.075
		spright=1-spleft; sptop=0.95; spbottom=0.275
		figsize = (8,11)
	else:
		hspace=0.05; wspace=0.05; spleft=0.05
		spright=0.95; sptop=0.95; spbottom=0.1
		figsize = (6,10)

	qlabelx=0.2; qlabely=0.045

	# Call plot maker
	pt.plotMbyN(COMPS, X=X, Y=Y, vmin=vmin, vmax=vmax, \
	u=UCOMPS, v=VCOMPS, Xd=Xd, Yd=Yd, \
	scale=scale, qs=qs, qlabelx=qlabelx, qlabely=qlabely, \
	spleft=spleft, spright=spright, sptop=sptop, spbottom=spbottom, \
	wspace=wspace, hspace=hspace, isf=iceC, \
	suptitle=suptitle, insetTitle=insetTitle, insetTitleX=insetTitleX, \
	cmap='coolwarm', contourfNlevels=contourfNlevels, \
	xticks=xticks, yticks=yticks, \
	xticklabels=xticklabels, yticklabels=yticklabels, \
	contourLevels=contourLevels, contour=contour, \
	contour2Levels=contour2Levels, contour2=contour2, \
	figsize=figsize, show=False, save=True, outname=f'Figure_4.png', \
	ts_data=ts_data, ts_time=ts_time, ts_lines=ts_lines, ts_labels=ts_labels, \
	ts_xlim=ts_xlim, ts_ylim=ts_ylim, ts_insetTitle=ts_insetTitle, \
	ts_legendLoc=3)

	quit()

#==

	

