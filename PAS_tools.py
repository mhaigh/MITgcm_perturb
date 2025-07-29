# PAS_tools.py

# tools.py

import sys

import matplotlib.pyplot as plt

from time import ctime

import numpy as np
from sklearn.linear_model import LinearRegression
from scipy.stats import pearsonr
from scipy.stats import t as tstudent
		
#====================================================================================

# Start/end month in PAS
TS_SHIFT = -2; TE_SHIFT = -3
T_START = 24*12 + 5*12 + TS_SHIFT
T_END = -11 + TE_SHIFT
IPO_T_START = 14*12+1 + TS_SHIFT
IPO_T_END = TE_SHIFT

# Fixed lons that correspond to lons used in computation of slope current etc.
WEST_GETZ_LONS = [1, 60]
GETZ_LONS = [60, 65]
WEST_PITW_LONS = [65, 105]
PITW_LONS = [105, 120]
WEST_PITE_LONS = [120, 180]
PITE_LONS = [180, 230]
EAST_PITE_LONS = [230, 250]
SLOPE_SECTIONS = {'westGetz':WEST_GETZ_LONS, 'Getz':GETZ_LONS, 'westPITW':WEST_PITW_LONS, 'PITW':PITW_LONS, 'westPITE':WEST_PITE_LONS, 'PITE':PITE_LONS, 'eastPITE':EAST_PITE_LONS}

#==

def getSectionIndices(sections):
	'''Return list of indices for given sections along continental slope.'''
	
	section_indices = []
	for section in sections:
		for xi in range(SLOPE_SECTIONS[section][0], SLOPE_SECTIONS[section][1]):
			section_indices.append(xi)
			
	return section_indices

#==

# General helper function to get the area of each cell from latitude and longitude arrays giving the coordinates of the cell centres. Adapted from Kaitlin Naughten's mitgcm_python.
def dA_from_latlon (lon, lat, periodic=False, return_edges=False):

	deg2rad = np.pi / 180.
	rEarth = 6.371e6

	# Make sure they're 2D
	if len(lon.shape) == 1 and len(lat.shape) == 1:
		lon, lat = np.meshgrid(lon, lat)
	# Now make the edges
	# Longitude
	if periodic:
		lon_extend = wrap_periodic(lon, is_lon=True)
		lon_edges = 0.5*(lon_extend[:,:-1] + lon_extend[:,1:])
	else:
		lon_edges_mid = 0.5*(lon[:,:-1] + lon[:,1:])
		# Extrapolate the longitude boundaries
		lon_edges_w = 2*lon_edges_mid[:,0] - lon_edges_mid[:,1]
		lon_edges_e = 2*lon_edges_mid[:,-1] - lon_edges_mid[:,-2]
		lon_edges = np.concatenate((lon_edges_w[:,None], lon_edges_mid, lon_edges_e[:,None]), axis=1)
	dlon = lon_edges[:,1:] - lon_edges[:,:-1] 
	# Latitude
	lat_edges_mid = 0.5*(lat[:-1,:] + lat[1:,:])
	lat_edges_s = 2*lat_edges_mid[0,:] - lat_edges_mid[1,:]
	lat_edges_n = 2*lat_edges_mid[-1,:] - lat_edges_mid[-2,:]
	lat_edges = np.concatenate((np.expand_dims(lat_edges_s,0), lat_edges_mid, np.expand_dims(lat_edges_n,0)), axis=0)
	dlat = lat_edges[1:,:] - lat_edges[:-1,:]
	# Now convert to Cartesian
	dx = rEarth*np.cos(lat*deg2rad)*dlon*deg2rad
	dy = rEarth*dlat*deg2rad
	dA = dx*dy
	if return_edges:
		return dA, lon_edges, lat_edges
	else:
		return dA
		
#==

def latlon_to_xy(X, Y, returnGridSpacing=True, deg2rad=np.pi/180., rEarth=6.371e6):
	'''Convert lat-lon (X,y) grid to Cartesian (x, y) grid using the Haversine formula.'''
	
	X *= deg2rad; Y *= deg2rad
	
	nY, nX = X.shape
	x = np.zeros((nY, nX)); y = np.zeros((nY, nX))
	dx = np.zeros((nY-1, nX-1)); dy = np.zeros((nY-1, nX-1))
	for i in range(1,nX):
		for j in range(1,nY):
			dX = X[j,i] - X[j,i-1]
			dY = Y[j,i] - Y[j-1,i]
			
			# Compute grid for x
			a = np.cos(Y[j,i])**2 * np.sin(dX/2)**2
			dx[j-1,i-1] = 2 * rEarth * np.arctan2(np.sqrt(a), np.sqrt(1-a)) 
			x[j,i] = x[j,i-1] + dx[j-1,i-1]
			
			# Repeat for y
			a = np.sin(dY/2)**2
			dy[j-1,i-1] = 2 * rEarth * np.arctan2(np.sqrt(a), np.sqrt(1-a)) 
			y[j,i] = y[j-1,i] + dy[j-1,i-1] 
	
	if returnGridSpacing:
		return x, y, dx, dy
	else:
		return x, y
		
#==

def getBearing_(lat1, lon1, lat2, lon2, deg2rad=np.pi/180.):
	
	lat1 = deg2rad * np.copy(lat1)
	lat2 = deg2rad * np.copy(lat2)
	lon1 = deg2rad * np.copy(lon1)
	lon2 = deg2rad * np.copy(lon2)

	dLon = (lon2 - lon1)

	y = np.sin(dLon) * np.cos(lat2)
	x = np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(dLon)

	brng = np.arctan2(y, x)

	brng = np.rad2deg(brng)

	return brng
	
#==

def getBearing(slope_x, slope_y):

	bearing = np.zeros(len(slope_x))
	bearing[0] = getBearing_(slope_y[0], slope_x[0], slope_y[1], slope_x[1])
	bearing[-1] = getBearing_(slope_y[-2], slope_x[-2], slope_y[-1], slope_x[-1])
	bearing[1:-1] = getBearing_(slope_y[:-2], slope_x[:-2], slope_y[2:], slope_x[2:])
	
	return bearing
	
#==

def getSlopeContour(bathy, X, Y, pathno, level=-1000):

	cs = plt.contour(X, Y, bathy, levels=[level])
	x = cs.collections[0].get_paths()[pathno].vertices[:,0][::-1]
	y = cs.collections[0].get_paths()[pathno].vertices[:,1][::-1]
	plt.clf()
	
	return x, y

#==

def getDecimalTime(t, t_start=None, t_end=None, PAS_correction=15):
	'''Given time array t of seconds since epoch, return decimalised time in units years.'''
	
	months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
	monthDec = {}
	
	for mi in range(len(months)):
		monthDec[months[mi]] = mi/12.
	
	year = [float(ctime(int(t_))[-4:])-PAS_correction + monthDec[ctime(int(t_))[4:7]] for t_ in t[t_start:t_end]]

	return year

#==

def getTimeYear(path, TIMESTART=1970):
	'''Load time data from PAS_time.npy in relevant run/path.
	Convert from seconds since TIMESTART to year in decimal format.'''
	
	try:
		t = np.load(path+'/PAS_time.npy')
	except:
		t = np.load('/home/michael/Documents/data/'+path+'/post/PAS_time.npy')
	#else:
	#	print('Error: PAS_tools.getTimeYear. No PAS_time.npy file at path.')
		
	year = np.array(getDecimalTime(t))
	year = year - year[0] + TIMESTART

	return year

#==

def demean(f, axis=0):
	'''Remove the mean from a field along the given axis.'''

	u = f.copy()
	
	up = np.zeros((u.shape))
	dim = up.shape[axis]
	
	mean = np.mean(u, axis=axis)

	if axis == 0:
		for i in range(dim):
			up[i, ] = u[i, ] - mean

	elif axis == 1:
		for i in range(dim):
			up[:, i, ] = u[:, i, ] - mean
	
	elif axis == 2:
		for i in range(dim):
			up[:, :, i] = u[:, :, i] - mean

	return up
	
#==

def desurf(f, axis=1):
	'''Subtract surface flow from all levels of field f.'''
	
	u = f.copy()
	
	up = np.zeros((u.shape))
	dim = up.shape[axis]

	if axis == 0:
		for i in range(dim):
			up[i, ] = u[i, ] - u[0, ]

	elif axis == 1:
		for i in range(dim):
			up[:, i, ] = u[:, i, ] - u[:, 0, ]
	
	elif axis == 2:
		for i in range(dim):
			up[:, :, i] = u[:, :, i] - u[:, :, 0]
	
	return up
	
#==
	
def corr(u, v, axis=0):
	'''Get the correlation between two fields along axis.'''

	# corr = cov(u) * cov(v) / (sd(u) * sd(v))

	ud = demean(u, axis)
	vd = demean(v, axis)

	return np.mean(ud * vd, axis=axis) / (np.mean(ud**2, axis=axis) * np.mean(vd**2, axis=axis))**0.5

#==

def corrMap(u, v):
	'''Return map of correlation between 3D field u(t,y,x) and v.
	v is either a timeseries or array with same dimensions as u.'''
	
	nT, nY, nX = u.shape
	
	corrOut = np.zeros((nY, nX))
	
	if len(v.shape) == 1:
		for yi in range(nY):
			for xi in range(nX):
				corrOut[yi,xi] = corr(u[:,yi,xi], v)
	else:
		for yi in range(nY):
			for xi in range(nX):
				corrOut[yi,xi] = corr(u[:,yi,xi], v[:,yi,xi])
		
	return corrOut
	
#==

def autocorr(u):
	'''Compute autocorrelation for timeseries u'''
	
	nT = len(u)
	
	r = np.ones(nT-2)
	for lag in range(1, nT-2):
		r[lag] = corr(u[lag:], u[:-lag])
	
	return r
	
#==

def eFolding(u):
	'''Compute e-folding timescale for timeseries u.'''
	
	if len(u.shape) > 1:
		print('Error: PAS_tools.eFolding. Input timeseries must be 1-dimensional.')
		quit()
	
	nT = len(u)

	lag = 1
	efold = 1
	while lag < nT-2:
		if corr(u[lag:], u[:-lag]) < 1. / np.e:
			efold = lag
			lag += nT
		lag += 1
		
	return efold

#==

def EffDoF(u):
	'''Return effective number of degrees of freedom for timeseries u.
	Use formulation of Leith [J. Appl. Meteor., 1973, p. 1066].'''
	
	efold = eFolding(u)
	
	dof = 0.5 * len(u) / efold

	return dof
	
#==

def pval_edof(r, data, predictor, mask=None): 
	'''After computing correlation value corr between data and predictor,
	recompute pval map accounting for effective degrees of freedom.'''
	
	sided = 2

	if mask is None:
		mask = True
		
	if mask:
		# Determine effective degrees of freedom
		dof = max(EffDoF(data), EffDoF(predictor))
		# Determine the value of the t-statistic 
		t = r * np.sqrt(dof) / np.sqrt(1-r**2)
		# Calculate the p-value using Student's-t cumulative distribution function
		pval = sided * tstudent.cdf(abs(t), dof)
	else:
		pval = 0

	return pval

#==

def pval_edof_2D(corr, data, predictor): 
	'''After computing correlation map corr between data and predictor,
	recompute pval map accounting for effective degrees of freedom.'''
	
	nT, nY, nX = data.shape
	
	efoldP = EffDoF(predictor)
	sided = 2

	pval = np.zeros((nY, nX))
	for j in range(nY):
		for i in range(nX):
			r = corr[j,i]
			# Determine effective degrees of freedom
			dof = max(EffDoF(data[:,j,i]), efoldP)
			# Determine the value of the t-statistic 
			t = r * np.sqrt(dof) / np.sqrt(1-r**2)
			# Calculate the p-value using Student's-t cumulative distribution function
			pval = sided * tstudent.cdf(abs(t), dof)
	# End loop
	
	return pval

#==

def corress(v1, v2, binwidth, nin):
	'''Adapted to Python from Paul Holland's Matlab code.
	Calculate the correlation coefficient (r) and give the significance value 
	(p) for two vectors, v1,v2, accounting for the effective sample size (n)
	due to autocorrelations r1,r2, following Bretherton et al. J. Climate 1999.
	
	function [r,2-p,n,r1,r2]=corress(v1,v2,binwidth,nin) 

	The input 'binwidth' allows for binning of input timeseries before the
	calculations. binwidth=1 gives no binning. Binning is one way to proceed for 
	timeseries that are highly autocorrelated, and hence oversampled, since the 
	Bretherton formula breaks down in that case. This is questionable, but in 
	practice the r,p,n are independent of bin width as long as n is not small.
	
	The input 'nin' allows the opportunity to specify the effective sample size, 
	thus over-riding the calculation. set nin=0 to calculate the sample size.'''
	
	# Check vectors have the same length.
	n1 = len(v1)
	n2 = len(v2)
	if n1 != n2:
		print('Error: PAS_tools.corress. Input vectors have different lengths.')
		quit()

	# Bin the timeseries to reduce oversampling
	nbins = int(np.ceil(n1/binwidth))
	v1bin = np.zeros(nbins)
	v2bin = np.zeros(nbins)
	binlo = 0
	binhi = binwidth
	for b in range(nbins):
		v1bin[b] = np.ma.mean(v1[binlo:binhi])
		v2bin[b] = np.ma.mean(v2[binlo:binhi])
		binlo = binlo + binwidth
		if binhi is not None:
			binhi = binhi+binwidth
			if binhi > n1:
				binhi = None

	# Calculate the lag-1 autocorrelation for both variables 
	r1 = corr(v1bin[:-1],v1bin[1:])
	r2 = corr(v2bin[:-1],v2bin[1:])

	# warning if negative autocorrelation detected
	if (r1<0) or (r2<0):
		print('Warning: negative autocorrelation detected')

	# warn user that if autocorrelation is large in either vector, the effective sample size formula below is invalid
	#if (max(r1,r2)>0.5)
	#  warning('large autocorrelation detected; formula invalid')
	#end

	if nin == -1:
		n = max(EffDoF(v1), EffDoF(v2)) + 2 # Plus 2 since 2 subtracted for other cases.
	elif nin == 0:
		# Use the autocorrelation to calculate the effective sample size 
		n = nbins * (1-abs(r1*r2)) / (1+abs(r1*r2))
	else:
		# Use specified effective sample size
		n = nin
	print(n)

	# calculate the correlation between the two variables 
	r = corr(v1bin, v2bin)

	# Determine the value of the t-statistic 
	t = r * np.sqrt(n-2) / np.sqrt(1-r**2)

	# Calculate the p-value using Student's-t cumulative distribution function
	sided = 2
	p = sided * tstudent.cdf(abs(t), n-2)
	
	return r, 2-p, n, r1, r2

#==

def pearson(u, v):
	'''Use scipy.pearson to compute the correlation.'''
	
	corr, pval = pearsonr(u, v)
	
	return corr

#==

def corr2(u, v, axis=0):

	uav = np.mean(u, axis=axis)
	vav = np.mean(v, axis=axis)
	
	return np.mean((u-uav) * (v-vav), axis=axis) / (np.mean((u-uav)**2, axis=axis) *  np.mean((v-vav)**2, axis=axis))**0.5
	
#==

def variance(u, axis=0):
	'''Variance of field u along axis.'''
	
	ud = demean(u, axis)
	
	return np.mean(ud**2, axis=axis)

#==

def cov(u, v, axis=0):
	'''Get the covariance between two fields along axis.'''

	ud = demean(u, axis)
	vd = demean(v, axis)

	return np.mean(ud * vd, axis=axis)
	
#==

def crossCorr(u, v, nl, nl2, P_VALUE=True):
	'''Get cross correlations for number of lags nl.'''
	
	if nl != 2 * nl2 + 1:
		print('PAS_tools.crossCorr error: nl and nl2 incompatible.')
		quit()
		
	corruv = np.zeros(nl)
	p_value = np.zeros(nl)
	
	#corruv[nl2] = corr(u,v)
	corruv[nl2], p_value[nl2] = pearsonr(u,v)
	
	for li in range(1,nl2+1):
		# Shift u ahead of v
		#corruv[nl2-li] = corr(u[li:], v[:-li])
		corruv[nl2-li], p_value[nl2-li] = pearsonr(u[li:], v[:-li])
		# Shift v ahead of u (u, v[li:])
		#corruv[nl2+li] = corr(u[:-li], v[li:])	
		corruv[nl2+li], p_value[nl2+li] = pearsonr(u[:-li], v[li:])	
	
	if P_VALUE:
		return corruv, p_value
		
	else:
		return corruv
	
#==

def movingAv(a, n=3):
	ret = np.ma.cumsum(a.copy(), dtype=float, axis=0)
	ret[n:] = ret[n:] - ret[:-n]
	return np.array(ret[n - 1:] / n)
	
#==

def windowAv(a, n, nBuffer=None, av=True):
	'''Returns running mean with window length of n of data array a.'''
	
	# Just return input array in this case.
	if n in [0, 1, None]:
		return a
		
	# Else compute running mean
	b = np.zeros(np.array(a).shape)

	if n % 2 == 0:
		n += 1

	if nBuffer == None:
		nBuffer = n
	else:
		nBuffer += 1 # For compatability with operations on n.
		
	if av:
		for i in range((n-1)//2, len(b)-(n-1)//2+1):
			#print(a[i-(n-1)//2:i+(n-1)//2].shape)
			b[i] = np.ma.mean(a[i-(n-1)//2:i+(n-1)//2], axis=0)
		return b[(nBuffer-1)//2:-(nBuffer-1)//2+1]
		
	else:
		return a[(nBuffer-1)//2:-(nBuffer-1)//2+1]
		
#==

def windowAvMonthly(f, nyears):
	'''Running-mean average over nyears, but averaging over each month individually.
	Parameter nyears should be odd.'''
	
	nT = f.shape[0]
	if nT % 12 != 0:
		print('Error PAS_tools.windowAvMonthly: input data f must represent integer number of years.')
		quit()
	
	nyearsTot = int(nT/12)
	
	f_out = f.copy()
	
	ny2 = int((nyears-1)/2)
	for year in range(ny2, nyearsTot-ny2):
		for mi in range(12):
			f_tmp = f[mi:][::12][year-ny2:year+ny2+1]
			f_out[12*year+mi] = np.mean(f_tmp, axis=0)

	nmonths = nyears * 12	
	return f_out[(nmonths-1)//2:-(nmonths-1)//2-1]
		
#==

def windowAvLoop(a, n):

	b = np.zeros(a.shape)
	nt, ny, nx = b.shape
	for j in range(ny):
		print(j)
		for i in range(nx):
			b[:,j,i] = windowAv(a[:,j,i], n)
	
	return b

#==

def windowCorr(u, v, window_lengths, t=None, return_uv=False):
	'''Compute correlations between time series u and v average applying 
	running means for a range of running-mean window lengths.'''

	uw = []; vw = []; tw = []
	corrw = np.zeros(len(window_lengths))
	print(len(u))
	
	if t is None:
		t = np.linspace(1, len(u), len(u))
		
	tw.append(t.copy()); uw.append(u); vw.append(v)
	corrw[0] = corr(u, v)
	
	for wi in range(1, len(window_lengths)):
		utmp = movingAv(u, n=wi)
		vtmp = movingAv(v, n=wi)
		
		corrw[wi] = corr(utmp, vtmp)
		
		uw.append(utmp); vw.append(vtmp)
		tw.append(movingAv(t, n=wi))
		
	if return_uv:
		return corrw, tw, uw, vw		
	else:
		return corrw
		
#==

def fft(u):
	'''Return fast Fourier transform of u and corresponding frequencies.'''
	
	Fu = np.fft.fftshift(np.fft.fft(u))
	freq = np.fft.fftshift(np.fft.fftfreq(len(u)))
	period = 1. / freq
	
	return Fu, freq, period

#==

def fftPower(Fu, nl2):
	'''From fft of u, get power for each absolute frequency.'''
	
	F = np.abs(np.copy(Fu))

	power = np.real(F[nl2:])
	power[1:] += np.real(F[:nl2][::-1])
	
	return power
	
#==

def crossCorrWindowAv(uin, vin, window_lengths, nl, nl2):
	'''Plot correlation between (u,v) for a range of lags AND running-mean window lengths.'''
	
	# For each running-mean window length, compute cross correlation for all lags.
	
	u = np.copy(uin); v = np.copy(vin)
	
	corrw = np.zeros((len(window_lengths), nl))
	p_value = np.zeros((len(window_lengths), nl))
	
	corrw[0], p_value[0] = crossCorr(u, v, nl, nl2)
	
	for wi in range(1, len(window_lengths)):
		utmp = movingAv(u, n=wi)
		vtmp = movingAv(v, n=wi)
		
		corrw[wi], p_value[wi] = crossCorr(utmp, vtmp, nl, nl2)
		
	return corrw, p_value

#==

def computeComposites(data, timeseries, thresholdScale=1.):
	'''Compute composites of data for timeseries.
	Variable thresholdScale defines scaling of standard deviation.
	Input data has three dimensions, with first representing time.
	Input timeseries has one time dimension.'''
	
	tseries = np.array(timeseries).copy()
	tseries = demean(tseries, axis=0)
	
	std = thresholdScale * (np.mean(tseries**2)**0.5)
	nT, nY, nX = data.shape
	
	compPos = np.zeros((nY, nX)); nPos = 0
	compNeg = np.zeros((nY, nX)); nNeg = 0
	print(np.mean(tseries))
	for ti in range(nT):
		if tseries[ti] > std:
			compPos	+= data[ti]
			nPos += 1
		elif tseries[ti] < -std:
			compNeg += data[ti]
			nNeg +=1
			
	compPos /= nPos
	compNeg /= nNeg
	
	return compPos, compNeg

#==

def computeJointComposites(data, timeseries, thresholdScale=1.):
	'''Compute joint composites of data for timeseries.
	Variable thresholdScale defines scaling of standard deviation.
	Input data has three dimensions, with first representing time.
	Input timeseries has one time dimension.'''
	
	timeseries1 = timeseries[0]
	timeseries2 = timeseries[1]
	
	try:
		nT1 = len(timeseries1)
		nT2 = len(timeseries2)
		if nT1 != nT2:
			print('PAS_tools.computeJointComposites error: two timeseries must have same length.')
	except:
		print('PAS_tools.computeJointComposites error: input variable timeseries must be list of two predictors with same length.')
		quit()
	
	std1 = thresholdScale * (np.mean(timeseries1**2)**0.5)
	std2 = thresholdScale * (np.mean(timeseries2**2)**0.5)
	
	nT, nY, nX = data.shape
	
	compPosPos = np.zeros((nY, nX)); nPosPos = 0
	compPosNeg = np.zeros((nY, nX)); nPosNeg = 0
	compNegPos = np.zeros((nY, nX)); nNegPos = 0
	compNegNeg = np.zeros((nY, nX)); nNegNeg = 0
		
	for ti in range(nT):
		if timeseries1[ti] > std1:
			if timeseries2[ti] > std2:
				compPosPos += data[ti]
				nPosPos += 1
			elif timeseries2[ti] < -std2:
				compPosNeg += data[ti]
				nPosNeg += 1
		elif timeseries1[ti] < -std1:
			if timeseries2[ti] > std2:
				compNegPos += data[ti]
				nNegPos += 1
			elif timeseries2[ti] < -std2:
				compNegNeg += data[ti]
				nNegNeg += 1
			
	print('Number of samples: ')
	print(nPosPos, nPosNeg, nNegPos, nNegNeg)
	compPosPos /= nPosPos
	compPosNeg /= nPosNeg
	compNegPos /= nNegPos
	compNegNeg /= nNegNeg
	
	return compPosPos, compPosNeg, compNegPos, compNegNeg
	
#==

def computeFourCompositeRanges(s1, e1, s2, e2, s3, e3, s4, e4, thresh=None, year=None, tseries=None, opt='decYear', plotCheck=False):
	'''Compute time ranges (returing time indices), 
	of four ranges for computing composites.
	Variable opt gives option to get indices from decimal year (decYear)
	or for searching for specific index in given range (range).
	decYear options requires decimal year timeseries.
	range option requires a std threshold and tseries data.'''
	
	# Get indices for four time periods from years
	if opt == 'decYear':
		if year is None:
			print('Error: PAS_tools.computeFourCompositeRanges')
			print('year data must be provided')
			quit()
		s1i = np.argmin(np.abs(s1-year))
		e1i = np.argmin(np.abs(e1-year))
		s2i = np.argmin(np.abs(s2-year))
		e2i = np.argmin(np.abs(e2-year))
		s3i = np.argmin(np.abs(s3-year))
		e3i = np.argmin(np.abs(e3-year))
		s4i = np.argmin(np.abs(s4-year))
		e4i = np.argmin(np.abs(e4-year))

	# Get indices for four time periods using std of timeseries.
	elif opt == 'range':
		if thresh is None or tseries is None:
			print('Error: PAS_tools.computeFourCompositeRanges')
			print('thresh and tseries data must be provided')
			quit()
		# Normalise timeseries.
		tmp = np.array(tseries).copy()
		tmp = demean(tmp, axis=0)
		tmp = tmp / np.mean(tmp**2)**0.5

		s1i = s1[0]+np.argmin(np.abs(tmp[s1[0]:s1[1]]-thresh))
		e1i = e1[0]+np.argmin(np.abs(tmp[e1[0]:e1[1]]-thresh))
		s2i = s2[0]+np.argmin(np.abs(tmp[s2[0]:s2[1]]+thresh))
		e2i = e2[0]+np.argmin(np.abs(tmp[e2[0]:e2[1]]+thresh))
		s3i = s3[0]+np.argmin(np.abs(tmp[s3[0]:s3[1]]-thresh))
		e3i = e3[0]+np.argmin(np.abs(tmp[e3[0]:e3[1]]-thresh))
		s4i = s4[0]+np.argmin(np.abs(tmp[s4[0]:s4[1]]+thresh))
		e4i = e4[0]+np.argmin(np.abs(tmp[e4[0]:e4[1]]+thresh))	

	if plotCheck:
		fs = 12
		std = np.mean((tmp-np.mean(tmp))**2)**0.5
		std = std * np.ones(tmp.shape)
		plt.figure(figsize=(12,5))
		plt.plot(year, tmp, color='k')
		plt.plot(year, std, color='k', linestyle='dashed')
		plt.plot(year, -std, color='k', linestyle='dashed')
		plt.ylim(-2,2)
		plt.xlim(min(year), max(year))
		ax = plt.gca()
		#ax.fill_between(year, std, tmp, where=tmp>std, color='indianred')	
		#ax.fill_between(year, tmp, -std, where=-std>tmp, color='cornflowerblue', alpha=0.8)
		ax.fill_between(year,-2,2,where=year<year[s1i], color='indianred', alpha=0.6)
		ax.fill_between(year,-2,2,where=(year>=year[e2i])&(year<year[s3i]), color='indianred', alpha=0.6)
		ax.fill_between(year,-2,2,where=(year>=year[s1i])&(year<year[e1i]), color='indianred')
		ax.fill_between(year,-2,2,where=(year>=year[s3i])&(year<year[e3i]), color='indianred')

		ax.fill_between(year,-2,2,where=(year<year[s2i])&(year>=year[e1i]), color='cornflowerblue', alpha=0.6)
		ax.fill_between(year,-2,2,where=(year>=year[e3i])&(year<year[s4i]), color='cornflowerblue', alpha=0.6)
		ax.fill_between(year,-2,2,where=(year>=year[s2i])&(year<year[e2i]), color='cornflowerblue')
		ax.fill_between(year,-2,2,where=(year>=year[s4i])&(year<year[e4i]), color='cornflowerblue')


		#plt.axvline(year[s1i],color='b');plt.axvline(year[e1i],color='r')
		#plt.axvline(year[s2i],color='b');plt.axvline(year[e2i],color='r')
		#plt.axvline(year[s3i],color='b');plt.axvline(year[e3i],color='r')
		#plt.axvline(year[s4i],color='b');plt.axvline(year[e4i],color='r')
		
		plt.grid();# plt.tight_layout();
		plt.title(r'Baroclinicity (units std. dev.) and composite ranges', fontsize=fs)	
		plt.xticks(fontsize=fs); plt.yticks(fontsize=fs)
		plt.savefig('fourCompsCheck')

	return s1i, e1i, s2i, e2i, s3i, e3i, s4i, e4i

#==

def seasonalData(data, year):
	'''Return dictionary of timeseries data averaged over four seasons.'''
	
	jani = np.argmin(year[:12]-np.floor(year[:12]))
	if (jani - 1) % 3 != 0:
		print('Adjust start so that data starts with 3-months of a season')
		quit()
	if len(year) % 12 != 0:
		print('Adjust end so that data spans integer number of years.')
		quit()

	nyears = int(len(year)/12)
	seasons = ['Summer', 'Autumn', 'Winter', 'Spring']
	shift = int((jani-1)/3)
	seasons = seasons[-shift:] + seasons[:-shift]
	
	if len(data.shape) > 1:
		# Assume data has x,y dims
		nt, ny, nx = data.shape
		season0 = np.zeros((nyears, ny, nx)); season1 = np.zeros((nyears, ny, nx))
		season2 = np.zeros((nyears, ny, nx)); season3 = np.zeros((nyears, ny, nx))
	else:
		season0 = np.zeros(nyears); season1 = np.zeros(nyears)
		season2 = np.zeros(nyears); season3 = np.zeros(nyears)
		
	for yi in range(nyears):
		season0[yi] = np.mean(data[12*yi+0:12*yi+3], axis=0)
		season1[yi] = np.mean(data[12*yi+3:12*yi+6], axis=0)
		season2[yi] = np.mean(data[12*yi+6:12*yi+9], axis=0)
		season3[yi] = np.mean(data[12*yi+9:12*yi+12], axis=0)
	
	data_seasonal = {seasons[0]:season0, seasons[1]:season1, seasons[2]:season2, seasons[3]:season3}
	
	return data_seasonal

#==

def seasonalDataIPO(data, year, IPO, thresholdScale=1.):
	'''Return dictionary of timeseries data averaged over four seasons.
	Also return seasonal timeseries averaged during positive and negative IPO.'''

	std = thresholdScale * (np.mean(IPO**2)**0.5)
			
	jani = np.argmin(year[:12]-np.floor(year[:12]))
	if (jani - 1) % 3 != 0:
		print('Adjust start so that data starts with 3-months of a season')
		quit()
	if len(year) % 12 != 0:
		print('Adjust end so that data spans integer number of years.')
		quit()

	nyears = int(len(year)/12)
	seasons = ['Summer', 'Autumn', 'Winter', 'Spring']
	shift = int((jani-1)/3)
	seasons = seasons[-shift:] + seasons[:-shift]
	
	# Assume data has x,y dims
	nt, ny, nx = data.shape
	seasonData = np.zeros((5, 3, ny, nx))
	nPos = np.zeros(4);	nNeg = np.zeros(4)
		
	last = 0
	for yi in range(nyears):
		for si in range(4):
			for mi in range(3):
				index = 12*yi + 3*si + mi
				seasonData[si, 0, ] += data[index, ]
				if IPO[index] > std:
					seasonData[si, 1, ] += data[index, ]
					nPos[si] += 1
				elif IPO[index] < - std:
					seasonData[si, 2, ] += data[index, ]
					nNeg[si] += 1
					
	seasonData[:,0] /= (0.25 * nt)

	seasonData[4,0] = np.mean(seasonData[:,0], axis=0)
	seasonData[4,1] = np.sum(seasonData[:,1], axis=0) / np.sum(nPos) - seasonData[4,0]
	seasonData[4,2] = np.sum(seasonData[:,2], axis=0) / np.sum(nNeg) - seasonData[4,0]
	
	for si in range(4):
		seasonData[si,1] /= nPos[si]
		seasonData[si,1] -= seasonData[si,0]
		#seasonData[si,1] -= np.mean(seasonData[:,0], axis=0)
		
		seasonData[si,2] /= nNeg[si]
		seasonData[si,2] -= seasonData[si,0]
		#seasonData[si,2] -= np.mean(seasonData[:,0], axis=0)
	
	return seasonData, seasons
	
#==

def monthlyDataComp(data, year, pred, forceMultiple12=False, DEMEAN=True, thresholdScale=1.):
	'''Return dictionary of timeseries data averaged over each month.
	Also return monthly timeseries averaged during positive and negative IPO.
	If DEMEAN, remove monthly mean from pos/neg IPO data so that they represent anomalies.'''

	print('BUG IN THIS CODE THAT DETERMINES FIRST MONTH TO ALWAYS BE JAN')
	
	mean = np.mean(pred)
	std = thresholdScale * (np.mean((pred-mean)**2)**0.5)
			
	jani = np.argmin(year[:12]-np.floor(year[:12]))
	if forceMultiple12:
		if len(year) % 12 != 0:
			print('Adjust end so that data spans integer number of years.')
			quit()

	nyears = int(len(year)/12)
	months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'July', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
	months = months[-int(jani):] + months[:-int(jani)]
	
	# Assume data has x,y dims
	nt, ny, nx = data.shape
	monthData = np.zeros((13, 3, ny, nx))
	nPos = np.zeros(12); nNeg = np.zeros(12)
		
	# For each year, for each month, add data to all-month sum.
	# Check if IPO pos/neg, and data to corresponding sum.
	for yi in range(nyears):
		for mi in range(12):
			monthData[mi, 0, ] += data[12*yi + mi, ]
			if (pred[12*yi + mi]-mean) > std:
				monthData[mi, 1, ] += data[12*yi + mi, ]
				nPos[mi] += 1
			elif (pred[12*yi + mi]-mean) < - std:
				monthData[mi, 2, ] += data[12*yi + mi, ]
				nNeg[mi] += 1
					
	print(nPos, nNeg)
	# Get mean by dividing by number of years.
	monthData[:,0] /= (nt / 12.)

	# Compute data averaged over all moonths.
	monthData[12,0] = np.mean(monthData[:,0], axis=0)
	monthData[12,1] = np.sum(monthData[:,1], axis=0) / np.sum(nPos)
	monthData[12,2] = np.sum(monthData[:,2], axis=0) / np.sum(nNeg)
	
	# Normalise monthly data, and remove monthly mean.
	for mi in range(12):
		monthData[mi,1] /= nPos[mi]
		monthData[mi,2] /= nNeg[mi]
	
	# Remove monthly mean from pos/neg IPO absolute values.
	if DEMEAN:
		monthData[:,1] -= monthData[:,0]
		monthData[:,2] -= monthData[:,0]
		
	return monthData, months

#==

def surfaceRadiation(qa, Ua, Ta, LW, SW, SSS, hi, hs, bathy, draft, nTiS=1000, T0=250, T1=270):
	'''Return all terms in surface radiation balance.
	Requires computation of surface temperature, also returned.'''
	
	qa = qa.copy()
	Ua = Ua.copy()
	Ta = Ta.copy()
	LW = LW.copy()
	SW = SW.copy()
	SSS = SSS.copy()
	hi = hi.copy()
	hs = hs.copy()
	
	aa1 = 2663.5
	aa2 = 12.537
	bb1 = 0.622
	bb2 = 1.0 - bb1
	Ppascals = 100000.0
	lnTEN = np.log(10)
	
	SEAICE_lhEvap = 2.5e6
	SEAICE_lhFusion = 3.34e5
	Ls = SEAICE_lhEvap + SEAICE_lhFusion # lhSublim
	
	rhoa = 1.2 # SEAICE_rhoAir
	rhos = 3.3e2 # SEAICE_rhoSnow
	ca = 1.005e3 # SEAICE_cpAir
	CDi = 1.75e-3 # SEAICE_dalton
	
	eps_s = 9.5e-1 # Sea-ice/snow emissivity. (Same val. for snow/ice.)
	sigma = 5.67e-8 # SEAICE_boltzmann
	
	SEAICE_dTempFrz_dS = -5.75e-2
	SEAICE_tempFrz0 = 9.01e-2
	celsius2K = 273.15

	Tf = SEAICE_dTempFrz_dS * SSS + SEAICE_tempFrz0 + celsius2K
	
	nY, nX = qa.shape
	TiS = np.zeros((nTiS, nY, nX))
	for yi in range(nY):
		for xi in range(nX):
			TiS[:,yi,xi] = np.linspace(T0, T1, nTiS)
	
	#==
	
	# LATENT
	
	mm_log10pi = -aa1 / TiS + aa2
	mm_pi = np.exp(mm_log10pi * lnTEN)
	qsat = bb1 * mm_pi / ( Ppascals - (1.0 - bb1) * mm_pi ) # qhice
	# qhice=qsat: saturation vapor pressure of snow/ice surface
	
	LATENT = rhoa * Ls * CDi * Ua * (qsat - qa)

	#==
	
	# SENS

	SENS = rhoa * ca * CDi * Ua * (TiS - Ta)

	#==
	
	# BLACKBODY
	
	BB = eps_s * sigma * TiS**4
	
	#==
	
	# LONGWAVE
	LWtmp = - eps_s * LW
	LW = np.zeros((nTiS, nY, nX))
	for ti in range(0,nTiS):
		LW[ti] = LWtmp
	
	#==
	
	# SHORTWAVE
	
	# DESIGNED TO WORK FOR WINTER ONLY, WHEN SEA ICE IS SNOW-COVERED.
	penetSWFrac = 0.0
	
	# These are values for DRY snow/ice, suitable for winter when Ta < surface melt temp.
	ALB_ICE = 7.5e-1
	ALB_SNOW = 8.4e-1
	
	# For albedo should compute a linear transition between ALB_SNOW and ALB_ICE where hsnow>HCUT=0.15.
	# But just using ALB_SNOW for now.
	ALB = ALB_SNOW
     
	SWtmp = - (1.0 - ALB) * (1.0 - penetSWFrac) * SW    
	SW = np.zeros((nTiS, nY, nX))
	for ti in range(0,nTiS):
		SW[ti] = SWtmp
		
	#==
	
	# CONDUCTIVE HEAT FLUX
	
	ki = 2.16560
	ks = 3.10000e-1
	
	FC = ki * ks * (Tf - TiS) / (ki * hs + ks * hi)
	
	#==
	
	# Find TiS that balances surface radiation:
	# Take argmin of poly and linearly interpolate.
	
	poly = LATENT + SENS + LW + BB + SW - FC
	
	SW_out = SW[0]; LW_out = LW[0]
	
	LATENT_out = np.zeros((nY, nX))
	SENS_out = np.zeros((nY, nX))
	BB_out = np.zeros((nY, nX))
	FC_out = np.zeros((nY, nX))
	
	tsurf = np.zeros((nY, nX))
	for yi in range(nY):	
		for xi in range(nX):
			p = poly[:, yi, xi]
			T = TiS[:, yi, xi]
			t_arg = np.argmin(np.abs(p))
			
			if bathy[yi,xi] < 0 and draft[yi,xi] == 0:
			
				if np.min(p)*np.max(p)>0:
					print('Polynomial has no solutions at (xi, yi) = ' + str((xi,yi)))
					plt.plot(T, p); plt.show()
									
				else:
					tsurf[yi,xi] = T[t_arg]
					LATENT_out[yi, xi] = LATENT[t_arg, yi, xi]
					SENS_out[yi, xi] = SENS[t_arg, yi, xi]
					BB_out[yi, xi] = BB[t_arg, yi, xi]
					FC_out[yi, xi] = FC[t_arg, yi, xi]

	# End loops
	#==
			
	return {'LATENT':LATENT_out, 'SENS':SENS_out, 'BB':BB_out, 'LW':LW_out, 'SW':SW_out, 'FC':FC_out, 'TiS':tsurf}

#==

def Fmi(SST, SSS):
	'''Compute flux from ocean mixed layer into base of sea ice.'''
	
	# Flux is 
	# f_mi = rhow * cw * ch * umix * (Tmix - Tf)
	
	rhow = 1026.
	cw = 4190.
	ch = 0.035 #0.006
	
	SEAICE_dTempFrz_dS = -5.75e-2
	SEAICE_tempFrz0 = 9.01e-2
	celsius2K = 273.15

	Tf = SEAICE_dTempFrz_dS * SSS + SEAICE_tempFrz0 + celsius2K
	Tmix = SST + celsius2K

	dzSurf = 10.
	deltaTtherm = 800.
	umix = dzSurf / deltaTtherm
	
	# https://github.com/MITgcm/MITgcm/blob/master/pkg/seaice/seaice_growth.F
	# https://github.com/MITgcm/MITgcm/blob/master/pkg/seaice/seaice_init_fixed.F
	# https://github.com/MITgcm/MITgcm/blob/master/pkg/seaice/SEAICE_PARAMS.h
	
	f_mi = rhow * cw * ch * umix * (Tmix - Tf)
	
	#SEAICE_frazilFrac*dzSurf/SEAICE_deltaTtherm
	#*
	#(HeatCapacity_Cp*rhoConst * recip_QI)
    # &         * (theta(i,j,kSurface,bi,bj)-tempFrz)
    # &         * SEAICE_deltaTtherm * HEFFM(i,j,bi,bj)
     
	return f_mi, rhow * cw * ch * umix * Tmix, -rhow * cw * ch * umix * Tf
	
#==

def surfaceRadiation_testCode():

	SW_out = SW[0]; LW_out = LW[0]
	
	LATENT_out = np.zeros((nY, nX))
	SENS_out = np.zeros((nY, nX))
	BB_out = np.zeros((nY, nX))
	FC_out = np.zeros((nY, nX))
	
	tsurf = np.zeros((nY, nX))
	tsurf2 = np.zeros((nY, nX))
	for yi in range(nY):	
		for xi in range(nX):
			p = poly[:, yi, xi]
			T = TiS[:, yi, xi]
			
			if bathy[yi,xi] < 0 and draft[yi,xi] == 0:
			
				if np.min(p)*np.max(p)>0:
					print(hs[yi,xi])
					print(hi[yi,xi])
					
					plt.subplot(121)
					plt.plot(T, p, label='poly')
					plt.plot(T, FC[:,yi,xi], label='FC')
					#plt.plot(T, SENS[:,yi,xi], label='SENS')
					#plt.plot(T, LW[:,yi,xi], label='LW')
					#plt.plot(T, SW[:,yi,xi], label='SW')
					#plt.plot(T, BB[:,yi,xi], label='BB')
					#plt.plot(T, LATENT[:,yi,xi], label='LATENT')
					plt.grid();	plt.legend()
					plt.subplot(122)
					plt.contourf(bathy, vmin=-1000, vmax=0, mesh=True); plt.plot(xi,yi, 'ro')
					plt.show()
									
				else:

					t_arg = np.argmin(np.abs(p))
					if p[t_arg] > 0:
						tsurf[yi,xi] = T[t_arg] + (T[t_arg] - T[t_arg-1]) / (p[t_arg] - p[t_arg-1])
					else:
						tsurf[yi,xi] = T[t_arg] - p[t_arg] * (T[t_arg] - T[t_arg+1]) / (p[t_arg] - p[t_arg+1])
					tsurf2[yi,xi] = T[t_arg]
					
					LATENT_out[yi, xi] = LATENT[t_arg, yi, xi]
					SENS_out[yi, xi] = SENS[t_arg, yi, xi]
					BB_out[yi, xi] = BB[t_arg, yi, xi]
					FC_out[yi, xi] = FC[t_arg, yi, xi]
			
	#pt.plot1by2([tsurf, tsurf2], vmin=[250]*2, vmax=[300]*2, mesh=True)
	
	RADS2 = {'LATENT':LATENT_out, 'SENS':SENS_out, 'BB':BB_out, 'LW':LW_out, 'SW':SW_out, 'FC':FC_out, 'TiS':tsurf2}

	return
	
#==

def monthlyData(data, year, forceMultiple12=False, correctMonthLag=False):
	'''Return mean of data for each month.'''
	
	print('BUG IN THIS CODE THAT DETERMINES FIRST MONTH TO ALWAYS BE JAN')
	nM = 12
	
	jani = np.argmin(year[:12]-np.floor(year[:12]))
	if forceMultiple12:
		if len(year) % 12 != 0:
			print('PAS_tools.monthlyData. Adjust end so that data spans integer number of years.')
			quit()
		
	nyears = int(len(year)/nM)
	months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'July', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
	if correctMonthLag:
		months = months[-1:] + months[:-1]
	months = months[-int(jani):] + months[:-int(jani)]
	
	if len(data.shape) == 1:
		nT = data.shape[0]
		dataMonthly = np.zeros((nM))
	elif len(data.shape) == 3:
		nT, nY, nX = data.shape
		dataMonthly = np.zeros((nM, nY, nX))

	if forceMultiple12:
		if nT % 12 != 0:
			print('Warning: PAS_tools.monthlyData. Length of first dim of data should be multiple of 12.')

	for mi in range(nM):
		dataMonthly[mi] = np.mean(data[mi:None:12], axis=0)
	
	return dataMonthly, months
	
#==

def monthlyTimeseries(data, year, forceMultiple12=False, correctMonthLag=False):
	'''Return mean of data for each month.'''
	
	print('BUG IN THIS CODE THAT DETERMINES FIRST MONTH TO ALWAYS BE JAN')
	nM = 12
	
	jani = np.argmin(year[:12]-np.floor(year[:12]))
	if forceMultiple12:
		if len(year) % 12 != 0:
			print('PAS_tools.monthlyData. Adjust end so that data spans integer number of years.')
			quit()
		
	nyears = int(len(year)/nM)
	months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'July', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
	if correctMonthLag:
		months = months[-1:] + months[:-1]
	months = months[-int(jani):] + months[:-int(jani)]
	
	if forceMultiple12:
		if nT % 12 != 0:
			print('Warning: PAS_tools.monthlyData. Length of first dim of data should be multiple of 12.')
	
	dataMonthly = []
	for mi in range(nM):
		print(len(data[mi:None:12]))
		dataMonthly.append(data[mi:None:12])
	
	return dataMonthly, months

#==

def deseason(uin):
	'''Deseasons monthly data. Checks if data has 12*n data points.'''
	
	u = np.copy(uin)
	
	nm = len(u)
	if nm % 12 != 0:
		print('Warning: PAS_tools.deseason. Length of uin should be multiple of 12.')
	
	cycle = np.zeros(12)
	for mi in range(12):
		cycle[mi] = np.mean(u[mi:None:12])
	
	cycle = np.tile(cycle, nm//12)
	
	if len(u.shape) > 1:
		return np.transpose(np.transpose(u, (2,1,0)) - cycle, (2,1,0))
	else:
		return u - cycle
		
	
#==

def detrendXY(uin, tin=None, interceptFlag=1):
	'''Use linear regression to detrend uin timeseries,
	where uin has dimensions of (time, lat, lon).'''	
	
	nt, ny, nx = uin.shape
	
	uout = uin.copy()
	
	if tin is None:
		tin = np.linspace(1,nt,nt)
	t = np.copy(tin).reshape(-1,1)
	
	for j in range(ny):
		for i in range(nx):
			reg = LinearRegression().fit(t, uin[:,j,i])
			uout[:,j,i] = uin[:,j,i] - (reg.coef_*tin + interceptFlag * reg.intercept_)
			
	return uout
	
#==

def detrendX(uin, tin=None, interceptFlag=1):
	'''Use linear regression to detrend uin timeseries,
	where uin has dimensions of (time, spatial dimension).'''	
	
	nt, nx = uin.shape
	
	uout = uin.copy()
	
	if tin is None:
		tin = np.linspace(1,nt,nt)
	t = np.copy(tin).reshape(-1,1)
	
	for i in range(nx):
		reg = LinearRegression().fit(t, uin[:,i])
		uout[:,i] = uin[:,i] - (reg.coef_*tin + interceptFlag * reg.intercept_)
			
	return uout
	
#==

def detrend(uin, tin=None, interceptFlag=1):
	'''Use linear regression to detrend uin timeseries.'''
	
	if len(uin.shape) == 2:
		return detrendX(uin, tin=tin, interceptFlag=interceptFlag)
	if len(uin.shape) == 3:
		return detrendXY(uin, tin=tin, interceptFlag=interceptFlag)
		
	# ELSE
	if tin is None:
		tin = np.linspace(1,len(uin),len(uin))
	
	u = np.copy(uin)
	t = np.copy(tin).reshape(-1,1)
	
	reg = LinearRegression().fit(t, u)
	
	return uin - (reg.coef_*tin + interceptFlag * reg.intercept_)
	
#==

def norm(data, printVarName=None, outputStats=False):
	'''Get time-mean and standard deviation of timeseries data.
	Return data normalised by mean and std.'''
	
	if len(data.shape) > 1:
		print('Error: PAS_tools.normaliseSTD. Input data must be 1D array.')
		
	dataAv = np.mean(data)
	dataSTD = np.mean((data-dataAv)**2)**0.5
	dataNorm = (data - dataAv) / dataSTD
	
	if printVarName is not None:
		print(printVarName)
		print('Mean: %.3f' % dataAv)
		print('STD: %.3f' % dataSTD)
	
	if outputStats:
		return dataNorm, dataAv, dataSTD
	else:
		return dataNorm
	
#==

def getCoriolis(lats, Omega=2*np.pi/86400., nx=None):
	'''Get Coriolis coefficient for spherical polar grid.
	Input lats to be in degrees.'''
	
	if len(lats.shape) == 2 or nx is None:
		f = 2. * Omega * np.sin(lats * np.pi / 180.)	
			
	elif len(lats.shape) == 1:
		f = np.zeros((lats.shape[0],nx))
		for i in range(nx):
			f[:,i] = 2. * Omega * np.sin(lats * np.pi / 180.)
	
	else:
		print('PAS_tools.getCoriolis error: lats must have 1 or 2 dimensions.')
		quit()
			
	return f

#==

def maskSouthNorthSlope(data, grid, lonsN=[225, 260], lonsS=[225,260], Nrange=58, FOR_PLOTTING=False):
	'''Return copies of data masked either side of continental slope.'''
	
	from plotting_tools import makeBathyContour
	
	bathyDeep = -1001
	bathyShallow = -999
	bathy = grid.bathy
	X = grid.XC; Y = grid.YC
	bathyC = makeBathyContour(bathy, X, Y)
	axis = 0
	
	if len(data.shape) > 2:
		nT, nY, nX = data.shape
		X = np.tile(X, (nT,1,1)).reshape((nT,nY,nX))
		Y = np.tile(Y, (nT,1,1)).reshape((nT,nY,nX))
		bathyC = np.tile(bathyC, (nT,1,1)).reshape((nT,nY,nX))
		axis = 1
		
	dataS = np.ma.masked_where(bathyC<bathyShallow, data)
	dataS = np.ma.masked_where(X<lonsS[0], dataS)
	dataS = np.ma.masked_where(X>lonsS[1], dataS)
		
	dataN = np.ma.masked_where(bathyC>bathyDeep, data) 
	dataN = np.ma.masked_where(np.roll(bathyC, Nrange, axis=axis)<bathyDeep, dataN)
	dataN = np.ma.masked_where(X<lonsN[0], dataN)
	dataN = np.ma.masked_where(X>lonsN[1], dataN)
	
	#==
	
	# This section edits masks, mainly for plotting purposes to show boundaries more clearly in contour plots.
	if FOR_PLOTTING:
	
		# Northern mask
		maskN = np.ma.getmask(dataN)
		for i in range(X.shape[1]):
			j = -1
			while j > -X.shape[0]:
				if maskN[j-1,i] and not maskN[j,i]:
					maskN[j,i] = True
					maskN[j+1,i] = True
					maskN[j+2,i] = True
					j -= X.shape[0]
				j -= 1
		dataN = np.ma.masked_where(maskN, dataN)
		
		# Southern mask
		maskS = np.ma.getmask(dataS)
		for i in range(X.shape[1]):
			j = -1
			while j > -X.shape[0]:
				if maskS[j,i] and not maskS[j-1,i]:
					maskS[j,i] = True
					maskS[j-1,i] = True
					maskS[j-2,i] = True
					#maskS[j-3,i] = True
					j -= X.shape[0]
				j -= 1
		dataS = np.ma.masked_where(maskS, dataS)
	
	#==
	
	return dataS, dataN	
	
#== 

def avSlopeSections(data, axis=1, sections=['westGetz', 'westPITW', 'westPITE']):
	'''Average data over slope sections.'''

	data_sum = np.zeros(data.shape[0])
	nx = 0
	for section in sections:
		iw = SLOPE_SECTIONS[section][0]
		ie = SLOPE_SECTIONS[section][1]
		data_sum += np.sum(data[:,iw:ie], axis=1)
		nx += ie-iw
	
	data_av = data_sum / nx
	
	return data_av
		
#==

def getSlopeWinds(path='', fname='slope_uw_av.npy', t_start=None, t_end=None, sections=['westGetz', 'westPITW', 'westPITE'], DETREND=True, AV=True, nn=60, interceptFlag=1):
	'''Load and process slope winds.'''

	u = np.load(path+fname)[t_start:t_end]
        
	if sections is not None:
		u = avSlopeSections(u, sections=sections)
        
	if DETREND:
		u = detrend(u, None, interceptFlag=interceptFlag)
        
	if AV:
		u = windowAv(u, n=nn)

	return u

#==

def getSurfCurrent(path='', fname='surf_uv_av.npy', t_start=None, t_end=None, sections=['westGetz', 'westPITW', 'westPITE'], DETREND=True, AV=True, nn=60, interceptFlag=1):
	'''Load surface current computed by PAS.py.'''
	
	u = np.load(path+fname)[t_start:t_end]
	
	if sections is not None:
		u = avSlopeSections(u, sections=sections)
	
	if DETREND:
		u = detrend(u, None, interceptFlag=interceptFlag)
	
	if AV:
		u = windowAv(u, n=nn)
		
	return u 
	
#==

def getUndercurrent(path='', fname='slope_uv_zmax.npy', t_start=None, t_end=None, sections=['westGetz', 'westPITW', 'westPITE'], DETREND=True, AV=True, nn=60, interceptFlag=1):#, t_start=T_START, t_end=T_END):
	'''Load undercurrent computed by PAS.py.'''
	
	u = np.load(path+fname)[t_start:t_end]
	
	if sections is not None:
		u = avSlopeSections(u, sections=sections)
	
	if DETREND:
		u = detrend(u, None, interceptFlag=interceptFlag)
	
	if AV:
		u = windowAv(u, n=nn)
		
	return u
	
#==

def getBarocl(path='', fname=['slope_uv_zmax.npy','surf_uv_av.npy'], t_start=None, t_end=None, sections=['westGetz', 'westPITW', 'westPITE'], DETREND=True, AV=True, nn=60, interceptFlag=1):
	'''Return difference between undercurrent and surface current computed by PAS.py.'''
	
	u = np.load(path+fname[0])[t_start:t_end] - np.load(path+fname[1])[t_start:t_end]
	
	if sections is not None:
		u = avSlopeSections(u, sections=sections)
	
	if DETREND:
		u = detrend(u, None, interceptFlag=interceptFlag)
	
	if AV:
		u = windowAv(u, n=nn)
		
	return u
	
#==

def coastalMask(grid, maskRad=5, draft=True, omitBurkeIsland=False):
	'''Return array of 1s for ocean with maskRad grid points of land.'''
	
	bathy = grid.bathy
	ny, nx = bathy.shape
	
	if omitBurkeIsland:
		X = grid.XC; Y = grid.YC
		mask = np.where(X>256, -1, np.where(X<254,-1,1))
		mask = np.where(Y>-72.5, -1, np.where(Y<-73.25,-1,mask))
		bathy = np.where(mask>0,-1000,bathy)
		
		mask = np.where(X>258, -1, np.where(X<256,-1,1))
		mask = np.where(Y>-74., -1, np.where(Y<-74.2,-1,mask))
		bathy = np.where(mask>0,-1000,bathy)
	
		mask = np.where(X>228, -1, np.where(X<227,-1,1))
		mask = np.where(Y>-74., -1, np.where(Y<-74.2,-1,mask))
		bathy = np.where(mask>0,-800,bathy)

	mask = np.zeros((ny, nx))
	for j in range(ny):
		for i in range(nx):
			if bathy[j,i] >= 0:
				mask[max(j-maskRad,0):min(j+maskRad,ny), max(i-maskRad,0):min(i+maskRad,nx)] = 1
			if bathy[j,i]<-1000:
				mask[j,i] = 0

	if draft: 
		for j in range(ny):
			for i in range(nx):
				if grid.draft[j,i] < 0:
					mask[max(j-maskRad,0):min(j+maskRad,ny), max(i-maskRad,0):min(i+maskRad,nx)] = 1
		mask = np.where(grid.draft<0, 0, mask)	
				
	mask = np.where(bathy>=0, 0, mask)
	mask[200:] = 0
	
	return mask
	
#==

def ASpolynyaMask(grid, maskRad=10):
	'''Return array of 1s in region of Amundsen Sea polynya, with zeros elsewhere.'''
	
	bathy = grid.bathy
	ny, nx = bathy.shape
	
	mask = np.zeros((ny, nx))
	for j in range(ny):
		for i in range(nx):
			if bathy[j,i] >= 0:
				mask[max(j-maskRad,0):min(j+maskRad,ny), max(i-maskRad,0):i] = 1
	
	X = grid.XC; Y = grid.YC
	mask = np.where(X>249.6, 0, np.where(X<248.2,0,mask))
	mask = np.where(Y>-73.1, 0, np.where(Y<-74.3,0,mask))
	
	return mask
	
#==

def sumOverIceShelves(isfw, grid, shelves=['PIG', 'THW', 'CRO', 'DOT', 'GETZe', 'COS', 'ABB', 'GETZw']):
	'''Sum ice-shelf melting over a number of ice shelves.
	Returns a 1-D timeseries in units kg/yr.'''

	secs_per_year = 86400. * 365.25
	
	areaIS = grid.RAC * (1-grid.hFacW[0]) * (1-grid.hFacS[0])
		
	isfwtot = np.zeros(isfw.shape[0])
	
	for shelf in shelves:
		melt = grid.getIceShelf(isfw.copy(), shelf)
		melt = secs_per_year * np.ma.sum(melt*areaIS, axis=(1,2))
		isfwtot += melt
		
	return isfwtot


