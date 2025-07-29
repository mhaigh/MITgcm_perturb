 # plotting_tools.py

import sys

import numpy as np

import matplotlib.pyplot as plt

from time import ctime, asctime, gmtime
from datetime import date
import dateutil.relativedelta

#==========================================================

def maskBathyAll(data, grid, color='grey', timeDep=False):
	'''Mask bathymetry everywhere in data field.
	If timeDep=True, assumes data is time-dependent with
	time on the first axis.'''

	# Make bathy and z 3D fields.
	bathy = grid.bathy
	z = grid.RC
	bathy = np.tile(bathy, (z.shape[0], 1, 1))
	z = np.tile(z, (1, bathy.shape[1], bathy.shape[2]))

	if timeDep:
		print('Note: creating 4D arrays for masking can be slow.')
		Nt = data.shape[0]
		bathy = np.tile(bathy, (Nt,1,1,1))
		z = np.tile(bathy, (Nt,1,1,1))

	return np.ma.array(data, mask=z<bathy)

#==

def maskBathyAllT(data, grid):
	'''Calls above function in temporal for loop, to get around high mem. demands.'''

	Nt = data.shape[0]
	for ti in range(Nt):
		data[ti] = maskBathyAll(data[ti], grid)

	return data

#==

def maskBathyXY(data, grid, zi, color='grey', subregion=False, lats=[], lons=[], timeDep=False, d=1):
	'''Function to be called before plotting to mask the bathymetry in (X,Y)-slice.'''
		
	# If masking in subregion, grid needs to account for it.
	if subregion == True:
		try:
			XC = grid.XC[lats[0]:lats[1]+1, lons[0]:lons[1]+1]
			YC = grid.YC[lats[0]:lats[1]+1, lons[0]:lons[1]+1]
			bathy = grid.bathy[lats[0]:lats[1]+1, lons[0]:lons[1]+1]
		except ValueError:
			print('Error: plotting_tools.maskBathyXY. If subregion set to True, both lons and lats need to be defined.')
			
	else:
		XC = grid.XC[::d,::d]
		YC = grid.YC[::d,::d]
		bathy = grid.bathy[::d,::d]
	
	# The depth of the slice.
	z = grid.RC.squeeze()[zi]
	# If data is time-dependent, assume time dimension is first and tile draft and z arrays.
	if timeDep:
		NT = data.shape[0]
		z = np.tile(z, (NT, 1, 1))
		bathy = np.tile(bathy, (NT, 1, 1))
		
	return np.ma.array(data, mask=z<bathy)
	
#==

def maskBathyXZ(data, grid, color='grey', yi=0, subregion=False, lons=[], depths=[], timeDep=False):
	'''Function to be called before plotting to mask the bathymetry in (X,Z)-slice.'''
	
	if (len(data.shape) != 2 and timeDep == False):
		print('Error: plotting_toolts.maskBathyXZ. If data is more than 2-dimensional, timeDep must be set to True.')
		sys.exit()
		
	# If masking in subregion, grid needs to account for it.
	if subregion == True:
		try:
			RC = grid.RC.squeeze()[depths[0]:depths[1]+1]
			XC = grid.XC[:, lons[0]:lons[1]+1]
			bathy = grid.bathy[:, lons[0]:lons[1]+1]
			
		except ValueError:
			print('Error: plotting_tools.maskBathyXZ. If subregion set to True, both lons and depths need to be defined.')
			
	else:
		RC = grid.RC.squeeze()
		XC = grid.XC
		bathy = grid.bathy
	
	# Make z a depth array with (y,z)-dimensions.
	z = np.tile(RC, (XC.shape[1], 1)).T
	
	# Make draft array with (y,z)-dimensions.
	bathy = np.tile(bathy[yi, ], (RC.shape[0], 1))
	
	# If data is time-dependent, assume time dimension is first and tile draft and z arrays.
	if timeDep:
		NT = data.shape[0]
		z = np.tile(z, (NT, 1, 1))
		bathy = np.tile(bathy, (NT, 1, 1))
		
	return np.ma.array(data, mask=z<bathy)
	
#==

def maskBathyYZ(data, grid, color='grey', xi=0, subregion=False, lats=[], depths=[], timeDep=False, sameIasNC=False):
	'''Function to be called before plotting to mask the bathymetry in (Y,Z)-slice.'''
	
	if (len(data.shape) != 2 and timeDep == False):
		print('Error: plotting_toolts.maskBathyYZ. If data is more than 2-dimensional, timeDep must be set to True.')
		sys.exit()
	
	# If masking in subregion, grid needs to account for it.
	if subregion == True:
		if sameIasNC:
			depths = [depths[0], depths[1]-1]
			lats = [lats[0], lats[1]-1]
		try:
			RC = grid.RC.squeeze()[depths[0]:depths[1]+1]
			YC = grid.YC[lats[0]:lats[1]+1,]
			bathy = grid.bathy[lats[0]:lats[1]+1,]
			
		except ValueError:
			print('Error: plottingtools.maskBathyYZ. If subregion set to True, lats and depths need to be defined.')
			
	else:
		RC = grid.RC.squeeze()
		YC = grid.YC
		bathy = grid.bathy
	
	# Make z a depth array with (y,z)-dimensions.
	z = np.tile(RC, (YC.shape[0], 1)).T
	
	# Make draft array with (y,z)-dimensions.
	bathy = np.tile(bathy[:,xi], (RC.shape[0], 1))

	# If data is time-dependent, assume time dimension is first and tile draft and z arrays.
	if timeDep:
		NT = data.shape[0]
		z = np.tile(z, (NT, 1, 1))
		bathy = np.tile(bathy, (NT, 1, 1))
		
	return np.ma.array(data, mask=z<bathy)

#==

def maskBathyYZ_allX(data, grid, timeDep=False):
	'''Calls maskBathyYZ for a range of x gridpoints.
	Makes the assumption that x is on the final axis.
	If timeDep this is axis=3, if not timeDep this is axis=2.'''

	data_masked = data.copy()

	if timeDep:
		Nx = data.shape[3]
	else:
		Nx = data.shape[2]

	for xi in range(Nx):
		data_masked[..., xi] = maskBathyYZ(data[..., xi], grid, xi=xi, timeDep=timeDep)

	return data_masked
	
#==

def maskDraftXY(data, grid, zi, color='grey', subregion=False, lats=[], lons=[], timeDep=False, inverse=False, d=1):
	'''Function to be called before plotting to mask the draft in (X,Y)-slice.
	If inverse, mask areas that are NOT ice shelf.'''
		
	# If masking in subregion, grid needs to account for it.
	if subregion == True:
		try:
			XC = grid.XC[lats[0]:lats[1]+1, lons[0]:lons[1]+1]
			YC = grid.YC[lats[0]:lats[1]+1, lons[0]:lons[1]+1]
			draft = grid.draft[lats[0]:lats[1]+1, lons[0]:lons[1]+1]
		except ValueError:
			print('Error: plotting_tools.maskBathyXY. If subregion set to True, both lons and lats need to be defined.')
			
	else:
		XC = grid.XC[::d,::d]
		YC = grid.YC[::d,::d]
		draft = grid.draft[::d,::d]
	
	# The depth of the slice.
	z = grid.RC.squeeze()[zi]

	# If data is time-dependent, assume time dimension is first and tile draft and z arrays.
	if timeDep:
		NT = data.shape[0]
		z = np.tile(z, (NT, 1, 1))
		draft = np.tile(draft, (NT, 1, 1))
		
	if inverse:
		return np.ma.array(data, mask=z<=draft)
	else:
		return np.ma.array(data, mask=z>draft)

#==

def maskDraftYZ(data, grid, color='grey', xi=10, subregion=False, lats=[], depths=[], timeDep=False):
	'''Function to be called before plotting to mask the ice shelf draft.
	If subregion is True, lats and depths lets grid know what the subregion is.'''
	
	if (len(data.shape) != 2 and timeDep == False):
		print('Error: plotting_toolts.maskDraftYZ. If data is more than 2-dimensional, timeDep must be set to True.')
		sys.exit()
	
	# If masking in subregion, grid needs to account for it.
	if subregion == True:
		try:
			RC = grid.RC.squeeze()[depths[0]:depths[1]+1]
			YC = grid.YC[lats[0]:lats[1]+1,]
			draft = grid.draft[lats[0]:lats[1]+1,]
			
		except:
			print('Error: plottingtools.maskDraftYZ. If subregion set to True, lats and depths need to be defined.')
			sys.exit()
			
	else:  
		RC = grid.RC.squeeze()
		YC = grid.YC
		draft = grid.draft
			
	# Make z a depth array with (y,z)-dimensions.
	z = np.tile(RC, (YC.shape[0], 1)).T
	
	# Make draft array with (y,z)-dimensions.
	draft = np.tile(draft[:,xi], (RC.shape[0], 1))
	
	# If data is time-dependent, assume time dimension is first and tile draft and z arrays.
	if timeDep:
		NT = data.shape[0]
		z = np.tile(z, (NT, 1, 1))
		draft = np.tile(draft, (NT, 1, 1))
	
	return np.ma.array(data, mask=z>draft)

#==

def maskBathyDraftSurf(data, grid, color='grey', subregion=False, lats=[], lons=[], timeDep=False):
	'''Combines the two functions that mask surface fields into a neater one.'''

	if len(data.shape) >= 3:
		timeDep = True
		
	data = maskBathyXY(data, grid, 0, color=color, subregion=subregion, lats=lats, lons=lons, timeDep=timeDep)
	data = maskDraftXY(data, grid, 0, color=color, subregion=subregion, lats=lats, lons=lons, timeDep=timeDep)
	
	return data

#==

def getTextData(TIME, t_format, xloc, yloc, color='k', short_ctime=True, ndays=1, lag=False, lagTime=0, PAS=False, overrideYear=1970):
	'''Get text data for animation plots.'''

	if t_format == 'month':
		tscale_val = 86400. * 30
		t = TIME / tscale_val
		text = [t_format + ' ' + str(int(tt)) for tt in t]

	elif t_format == 'day':
		tscale_val = 86400. * ndays
		t = TIME / tscale_val
		if lag:
			text = ['lag = ' + str(int(lagTime-ndays*tt)) + ' days' for tt in t]
		else:
			text = [t_format + ' ' + str(int(ndays*tt)) for tt in t]
	
	elif t_format == 'nc':
		from netCDF4 import num2date
		dates = num2date(TIME,'seconds since 1971-01-01 00:00:00', 'Gregorian')
		if isinstance(dates, list):
			text = [date.strftime()[:10] for date in dates]
		else:
			text = [dates.strftime()[:10]]

	elif t_format == 'ctime' or t_format == 'gmtime' or t_format == 'datetime':
		if t_format == 'ctime':
			text = [ctime(int(TIME[ti])) for ti in range(len(TIME))]
		elif t_format == 'gmtime':
			text = [asctime(gmtime(int(TIME[ti]))) for ti in range(len(TIME))]
		else:
			text = []
			for time_ in TIME:
				d = date.fromtimestamp(int(time_))
				d = d - dateutil.relativedelta.relativedelta(months=1)
				text.append(d.ctime()) 
			
		if PAS:
			correction = -15
		else:
			correction = overrideYear - 1970

		if short_ctime:
			text = [t[4:7] + ' ' + str(int(t[-4:])+int(correction)) for t in text]
			print(text)
			
	elif t_format == 'X':
		text = ['x = ' + str(time_/1000) + ' km' for time_ in TIME]

	text_data = {'text':text, 'xloc':xloc, 'yloc':yloc, 'fontdict':{'fontsize':14, 'color':color}}
	
	return text_data

#==

def setText(ax, text_data, i=None, set_invisible=False):
	'''Utility function for setting text on axis given text_data dictionary.'''

	if set_invisible:
		for t in ax.texts:
			t.set_visible(False)
	
	if i is not None:
		ax.text(text_data['xloc'], text_data['yloc'], text_data['text'][i], fontdict=text_data['fontdict'], zorder=14)	
		
	else:
	
		if isinstance(text_data['xloc'], list):
			for ti in range(len(text_data['text'])):
				ax.text(text_data['xloc'][ti], text_data['yloc'][ti], text_data['text'][ti], fontdict=text_data['fontdict'][ti], zorder=14)
		else:
			ax.text(text_data['xloc'], text_data['yloc'], text_data['text'], fontdict=text_data['fontdict'], zorder=14)	

#==

def getContourfLevels(vmin, vmax, contourfNlevels):
	'''Get kwarg levels for contourf plot from given vmin and vmax.
	vmin and vmax can be None, in which case None is returned.'''

	if vmin == None or vmax == None:
		return None
	
	else:
		return np.linspace(vmin, vmax, contourfNlevels)
		
#==

def makeList(arg, m, n=None, FORCE_MAKE_LIST=False):
	'''Check if argument is list.
	If not list, make it a list of repeated argument.'''

	if not isinstance(arg, list) or FORCE_MAKE_LIST:
		if n is None:
			arg = [arg]*m
		else:
			arg = [[arg]*m]*n
			
	return arg

#==

def doTicks(xticks, xticksvis, yticks, yticksvis, fontsize=None):

	if xticks is not None:
		if xticksvis:				
			plt.xticks(xticks, fontsize=fontsize)
		else:			
			plt.xticks(xticks, labels='', fontsize=fontsize)
	if yticks is not None:
		if yticksvis:				
			plt.yticks(yticks, fontsize=fontsize)
		else:			
			plt.yticks(yticks, labels='', fontsize=fontsize)

#==

def doXticks(xticks, xticksvis, fontsize=None):

	if xticks is not None:
		if xticksvis:				
			plt.xticks(xticks, fontsize=fontsize)
		else:			
			plt.xticks(xticks, labels='', fontsize=fontsize)
			
#==

def doYticks(yticks, yticksvis, fontsize=None):
	
	if yticks is not None:
		if yticksvis:				
			plt.yticks(yticks, fontsize=fontsize)
		else:			
			plt.yticks(yticks, labels='', fontsize=fontsize)
			
#==

def doLabels(xlabel, ylabel, fontsize=12):

	if xlabel is not None:
		plt.xlabel(xlabel, fontsize=fontsize)
	if ylabel is not None:
		plt.ylabel(ylabel, fontsize=fontsize)

#==

def doXlabels(xlabel, fontsize=12):

	if xlabel is not None:
		plt.xlabel(xlabel, fontsize=fontsize)

#==

def doYlabels(ylabel, fontsize=12):

	if ylabel is not None:
		plt.ylabel(ylabel, fontsize=fontsize)
	
#==

def doTitle(title, fontsize=12):

	if title is not None:
		plt.title(title, fontsize=fontsize)
		
#== 

def doStippling(data, X=None, Y=None, stipData=[None,1,1,1], marker='o', color='k', stipLim2=None):

	dataLim = stipData[0]
	dx = stipData[1]
	dy = stipData[2]
	s = stipData[3]

	ny, nx = data.shape
	if X is None:
		X = np.linspace(0, nx, nx)
	if Y is None:
		Y = np.linspace(0, ny, ny)
	
	data = data[::dy, ::dx]
	
	if len(X.shape) > 1:
		X = X[0, ::dx]
	else:
		X = X[::dx]
		
	if len(Y.shape) > 1:
		Y = Y[::dy, 0]
	else:
		Y = Y[::dy]
		
	if stipLim2 is not None:
		ys, xs = np.where(data<=stipLim2) #* np.where(data>dataLim)
		ys = Y[ys]; xs = X[xs]
		plt.scatter(xs, ys, s=s, color='grey', marker=marker)
		
	#if dataLim is not None:
	#	data = np.where(data<=dataLim, 1, 0) 
	ys, xs = np.where(data<=dataLim)
	ys = Y[ys]; xs = X[xs]
	
	plt.scatter(xs, ys, s=s, color=color, marker=marker)

		
	# End
	
#==

def doBox(box, X, Y):
	'''If box is not None, assumes it's a list/tuple defining the edges of box to plot over data.
	Format is (ys, yn, xw, xe).'''

	if box is not None:
		
		ys = box[0]; yn = box[1]
		xw = box[2]; xe = box[3]
		
		ymax = np.max(Y); xmax = np.max(X)
		ymin = np.min(Y); xmin = np.min(X)
		
		ysl = (ys - ymin) / (ymax - ymin)
		ynl = (yn - ymin) / (ymax - ymin)
		xwl = (xw - xmin) / (xmax - xmin)
		xel = (xe - xmin) / (xmax - xmin)
		
		plt.axvline(x=xw, ymin=ysl, ymax=ynl, color='k', linewidth=1.0, linestyle='--')
		plt.axvline(x=xe, ymin=ysl, ymax=ynl, color='k', linewidth=1.0, linestyle='--')
			
		plt.axhline(y=ys, xmin=xwl, xmax=xel, color='k', linewidth=1.0, linestyle='--')
		plt.axhline(y=yn, xmin=xwl, xmax=xel, color='k', linewidth=1.0, linestyle='--')
		
#==

def makeBathyContour(bathy, X, Y):
	'''Crude and ad hoc method for isolating -1000m slope contour.
	For use in plotting slope contour.'''
	
	bathyC = bathy.copy()
	
	bathyC = np.where(Y<-74.1, -400, bathyC)
	bathyC = np.where(Y>-70., -2000, bathyC)

	#bathyC = np.where(X<225., bathy, bathyC)
	bathyC = np.where(X>270., bathy, bathyC)
		
	xw = 230
	yn = -73
	yni = int(np.argmin(np.abs(Y[:,0]-yn)))
	xni = int(np.argmin(np.abs(X[yni,:]-xw)))
	bathyC[:yni, xni:] = -400
	
	xw = 260
	yn = -72
	yni = int(np.argmin(np.abs(Y[:,0]-yn)))
	xni = int(np.argmin(np.abs(X[yni,:]-xw)))
	bathyC[:yni, xni:] = -400
		
	return bathyC
	
#==
	
def doAntarcticInset(fig, data):
#def doAntarcticInset(ax, data):
	'''Plot inset with map of Antarctica.'''

	from mpl_toolkits.basemap import Basemap

	#[0.104, 0.61, 0.2, 0.2]
	#[0.11, 0.66, 0.15, 0.15]
	#[0.0865, 0.651, 0.175, 0.175]
	#[0.107, 0.781, 0.175, 0.175]
	left, bottom, width, height = [0.08, 0.799, 0.16, 0.16]
	ax2 = fig.add_axes([left, bottom, width, height])
	
	#iax = ax.inset_axes([0, .8, .3, .3])
	#iax.set_aspect('equal', anchor="NW")
	
	# (latitude of true scale) is pole.
	m = Basemap(projection='spstere',boundinglat=-61.5,lon_0=180,resolution='l')#, ax=iax, anchor='NW')
	m.drawcoastlines()
	m.fillcontinents(color='darkgrey',lake_color='w')
	m.drawmapboundary(fill_color='aqua')
	NP = 15
	c1 = 'k'; c2 = 'grey'; lw = 0.7

	xw = -140; xe = -80; ys = -75.5; yn =  -62
	m.plot(np.linspace(xw, xe, NP), np.ones(NP)*ys, color=c1, linewidth=lw, latlon=True)
	m.plot(np.linspace(xw, xe, NP), np.ones(NP)*yn, color=c1, linewidth=lw, latlon=True)
	m.plot(np.ones(NP)*xw, np.linspace(ys,yn,NP), color=c1, linewidth=lw, latlon=True)
	m.plot(np.ones(NP)*xe, np.linspace(ys,yn,NP), color=c1, linewidth=lw, latlon=True)
	xw = 225; xe = 262; ys = -75.5; yn = -68
	m.plot(np.linspace(xw, xe, NP), np.ones(NP)*ys, color=c2, linewidth=lw, latlon=True)
	m.plot(np.linspace(xw, xe, NP), np.ones(NP)*yn, color=c2, linewidth=lw, latlon=True)
	m.plot(np.ones(NP)*xw, np.linspace(ys,yn,NP), color=c2, linewidth=lw, latlon=True)
	m.plot(np.ones(NP)*xe, np.linspace(ys,yn,NP), color=c2, linewidth=lw, latlon=True)
	
#==

def makeStaggeredGrid(u, dx, dy):
	'''Return u with data staggered at each latitude, for better quiver plots.'''
	
	print('Calling makeStaggeredGrid. Edits may be needed if dy > 2.')

	nY, nX = u[::dy,::dx].shape
	nY -= 1; nX -=1
	ud = np.zeros((nY,nX))
	for j in range(nY):
		if j % 2 == 0:
			ud[j,:] = u[dy*j,0:-1:dx]
		else:
			ud[j,:] = u[dy*j,int(dx/2):None:dx] 
			
	return ud
	
#==
