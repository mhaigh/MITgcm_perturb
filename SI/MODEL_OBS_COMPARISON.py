import scipy.io as spio
import scipy.interpolate as interp

import numpy as np
import matplotlib.pyplot as plt

from datetime import datetime as dt

from grid_PAS import Grid as Grid_PAS
import plotting as pt
import PAS_tools as ptools
import tools as ts

#=================================================

def loadmat(filename):
	'''
	this function should be called instead of direct spio.loadmat
	as it cures the problem of not properly recovering python dictionaries
	from mat files. It calls the function check keys to cure all entries
	which are still mat-objects
	'''
	data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True)
	return _check_keys(data)

def _check_keys(dict):
	'''
	checks if entries in dictionary are mat-objects. If yes
	todict is called to change them to nested dictionaries
	'''
	for key in dict:
		if isinstance(dict[key], spio.matlab.mio5_params.mat_struct):
			dict[key] = _todict(dict[key])
	return dict        

def _todict(matobj):
	'''
	A recursive function which constructs from matobjects nested dictionaries
	'''
	dict = {}
	for strg in matobj._fieldnames:
		elem = matobj.__dict__[strg]
		if isinstance(elem, spio.matlab.mio5_params.mat_struct):
			dict[strg] = _todict(elem)
		else:
			dict[strg] = elem
	return dict
	
#=================================================

# Data required:
# 1. PIGmooringHOVMOLLER_v2
# 2. THETA_X258_Y-74p8.npy
# 3. PAS_time.npy

path = 'data/SI/'
data = loadmat(path+'PIGmooringHOVMOLLER_v2.mat')['PIGS_grid']
# dtype=[('depth', 'O'), ('date', 'O'), ('temp', 'O'), ('Tfit', 'O'), ('Tcleanfit', 'O'), ('Scleanfit', 'O'), ('PTcleanfit', 'O')])}

time = np.array(data['date'], dtype=float) - (1970*365.25-15)
Z = -np.array(data['depth'], dtype=float)
temp = np.array(data['Tcleanfit'], dtype=float)
nZ, nT = temp.shape

time = 86400*time
# Lnes for testing times
#fts = dt.fromtimestamp
#print(time_tmp)
#print(time[0], time[0]*86400)
#print(0, fts(0))
#print(time_tmp[0], fts(time_tmp[0])); quit()
for ti in range(len(time)):
	time[ti] = ptools.toYearFraction(time[ti])
	
GET_MONTHLY = True
if GET_MONTHLY:
	
	# Get monthly data with two methods
	
	# First method, monthly running mean.
	time = ptools.windowAv(time, n=30)
	temp = ptools.windowAv(temp.T, n=30).T
	nZ, nT = temp.shape
	
	# Second method, sum days in each month.
	# 9 days in first January (2009).
	days = [9,28,31,30,31,30,31,31,30,31,30,31]
	daysNormalYear = [31,28,31,30,31,30,31,31,30,31,30,31]
	daysLeapYear = [31,29,31,30,31,30,31,31,30,31,30,31]
	fourYears = daysLeapYear+daysNormalYear+daysNormalYear+daysNormalYear
	# Add 2010, 2011.
	days = days+daysNormalYear+daysNormalYear
	days = days+fourYears+fourYears+fourYears+fourYears+fourYears+fourYears
	
	ti = 0
	nMonths = len(days)
	tempMonthly = np.zeros((nZ, nMonths))
	timeMonthly = []
	for mi in range(nMonths):
		ndays = days[mi]
		if ti+ndays < nT:
			tempMonthly[:,mi] = temp[:,ti:ti+ndays].mean(axis=1)
			ti += ndays
	tempMonthly = np.array(tempMonthly)
	timeMonthly = np.array([2009+mi/12 for mi in range(nMonths)])
	t_end = np.argmin(np.abs(timeMonthly-2024))
	# Make it run to 2024.
	tempMonthly = tempMonthly[:,:t_end]
	timeMonthly = timeMonthly[:t_end]
	nTmonthly = len(timeMonthly)
			
#dt = time_tmp[1] - time_tmp[0]
#nn = int(86400.*30 / dt)

#==

# PAS output
pathr = 'data/grid/'
pathp = 'data/'

# Grid
grid = Grid_PAS(pathr)
X = grid.XC; Y=grid.YC; Z_PAS = grid.RC.squeeze()
draft = grid.draft
bathy = grid.bathy
bathy = np.ma.masked_where(bathy>=0, bathy)
#bathy = np.ma.masked_where(draft<0, bathy)

#pt.plot1by1(bathy, vmin=-1100, vmax=-800, mesh=True)#, X=X, Y=Y,)

time_PAS = np.load(pathp+'PAS_time.npy')
nT_PAS = len(time_PAS)
for ti in range(nT_PAS):
	time_PAS[ti] = ptools.toYearFraction(int(time_PAS[ti]))+9-1./12
	
# Data
fname = 'THETA_X258_Y-74p8.npy'
#fname = 'THETA_X258_Y-75p2.npy'
temp_PAS = np.load(pathp+fname)

#==

# Get isotherm heights

THERM = 0

# Convective periods to be highlighted on plot
isotherm_PAS = np.zeros(nT_PAS)
for ti in range(nT_PAS):
	isotherm_PAS[ti] = ts.getIsothermHeight1D(temp_PAS[ti], THERM, Z_PAS, upperMask=-100)
	#if np.max(temp_PAS[ti]) <= 0:
	#	isotherm_PAS[ti] = np.nan
	tt = time_PAS[ti]
	if (2013.8<tt and tt<2014.3) or (2015.7<tt and tt<2016.3) or (1998.8<tt and tt<1999.1):
		isotherm_PAS[ti] = np.nan
		
#isotherm_PAS = np.ma.masked_where(np.max(temp_PAS,axis=1)<0, isotherm_PAS)
isotherm = np.zeros(nT)
for ti in range(nT):
	isotherm[ti] = ts.getIsothermHeight1D(temp[:,ti], THERM, Z, upperMask=-100)

isothermMonthly = np.zeros(nTmonthly)
for ti in range(nTmonthly):
	isothermMonthly[ti] = ts.getIsothermHeight1D(tempMonthly[:,ti], THERM, Z, upperMask=-100)
	
#==

# Correlation between isotherm heights

# Mooring time must be completely inside PAS time.
time_PAS_max = max(time_PAS)
t_end = np.argmin(np.abs(time-time_PAS_max))-1 

# Interpolate modelled isotherm onto mooring time array.
f = interp.interp1d(time_PAS, isotherm_PAS)
isotherm_PAS_interp = f(time[:t_end])

# Make two timeseries without nans for correlations
ts1 = []; ts2 = []
for ti in range(len(isotherm_PAS_interp)):
	if isotherm_PAS_interp[ti] == isotherm_PAS_interp[ti]:
		ts1.append(isotherm_PAS_interp[ti])
		ts2.append(isotherm[ti])

binwidth = 1
print(ptools.corress(np.array(ts1), np.array(ts2), binwidth, -1))
# 0.4055622423706731, 0.08451763006498769

# Repeat correlation for properly done monthly data.

# Mooring time must be completely inside PAS time.
t_start = np.argmin(np.abs(time_PAS-2009))
t_end = np.argmin(np.abs(timeMonthly-time_PAS_max))+1

#print(len(timeMonthly[:t_end]), len(time_PAS[t_start:]))
#plt.plot(time_PAS[t_start:]+0.1)
#plt.plot(timeMonthly[:t_end])
#plt.plot(); plt.show(); quit()

# Make two timeseries without nans for correlations
ts1 = []; ts2 = []
ts1_tmp = isotherm_PAS[t_start:]
ts2_tmp = isothermMonthly[:t_end]

for ti in range(len(ts1_tmp)):
	if ts1_tmp[ti] == ts1_tmp[ti]:
		ts1.append(ts1_tmp[ti])
		ts2.append(ts2_tmp[ti])
		
print(ptools.corress(np.array(ts1), np.array(ts2), binwidth, -1))
# 0.4055622423706731, 0.08451763006498769
	
#==

# Plot

vmin = -2; vmax = 1.5

temp_PAS = ts.boundData(temp_PAS, vmin, vmax)
temp = ts.boundData(temp, vmin, vmax)

tmin = max(time[0], time_PAS[0])
tmax = min(time[-1], time_PAS[-1]) 
xlim = (tmin, tmax)
ylim = (-670, -340)

plt.figure(figsize=(8,4))
plt.plot(time_PAS, isotherm_PAS, label='Model', color='k')
#plt.plot(time[:t_end], isotherm_PAS_interp, label='Model interp')
#plt.plot(timeMonthly, isothermMonthly, label='Mooring', color='red')
plt.plot(time, isotherm, label='Mooring', color='cornflowerblue')
plt.xlim(min(time_PAS), max(time))
plt.ylim(ylim)
plt.text(1980, -635, 'r=0.45, p=0.08', fontsize=14)

# CTDs
lw = 2.5; marker = 'gx--'
plt.plot(1994.0, -401.11, marker, linewidth=2, markersize=12)
#plt.plot(2007.1527, -370.3, 'ro--', linewidth=2, markersize=12)
plt.plot(2009.059874, -371.7750974, marker, linewidth=lw, markersize=12)
plt.plot(2010.16172, -408.0564491, marker, linewidth=lw, markersize=12)
plt.plot(2012.148537, -563.4932852, marker, linewidth=lw, markersize=12)
plt.plot(2014.120738, -484.020553, marker, linewidth=lw, markersize=12)
plt.plot(2016.089625, -503.3291661, marker, linewidth=lw, markersize=12)
plt.plot(2017.125347, -522.7106024, marker, linewidth=lw, markersize=12)
plt.plot(2019.180746, -419.2919172, marker, linewidth=lw, markersize=12)
plt.plot(2020.108191, -422.4635305, marker, linewidth=lw, markersize=12)
# End CTDs


# Highlight convective periods
where = isotherm_PAS!=isotherm_PAS
ti = 1
while ti < len(where)-1:
	if not where[ti-1]:
		if where[ti]:
			where[ti] = False
			ti += 1
	elif (where[ti] and not where[ti+1]):
		where[ti] = False
		ti += 1
	ti += 1
			
ax = plt.gca()
ax.fill_between(time_PAS, ylim[0], ylim[1], where=where, alpha=0.2, color='firebrick')


plt.legend()
plt.grid()
plt.title('PIG south 0-deg. isotherm depth (m)', fontsize=16)
plt.savefig('Figure_SI_mooring.png')
plt.show()

quit()

plt.figure(figsize=(7,4))
plt.contourf(time_PAS, Z_PAS, temp_PAS.T, vmin=vmin,vmax=vmax, cmap='coolwarm')
plt.xlim(xlim); plt.ylim(ylim)
plt.colorbar()
plt.plot(time_PAS, isotherm_PAS, color='yellow')
plt.title('Modelled PIG south pot. temp. (deg.C)')
plt.show()

plt.figure(figsize=(7,4))
plt.contourf(time, Z, temp, vmin=vmin,vmax=vmax, cmap='coolwarm')
plt.xlim(xlim); plt.ylim(ylim)
plt.colorbar()
plt.plot(time, isotherm, color='yellow')
plt.title('Observed PIG south cons. temp. (deg.C)')
plt.show()
quit()

# CTD data and 0-deg isotherm depths (in order): 
# depths = [-401.11, -460.83, -489.1, -460.87, -489.03, 
# -349.62, -348.22,-356.4,-355.95, -381.86, -380.89, -370.47, -370.3]

#1994.1985901826483 -75.061 -101.848
#1994.1991342909691 -74.904 -101.474
#1994.199505327245 -74.741 -101.365
#1994.1991342909691 -74.904 -101.474
#1994.199505327245 -74.741 -101.365
#2007.1527111872147 -74.8005 -102.9453
#2007.1528234081684 -74.80059814453125 -102.94539642333984
#2007.1530954781836 -74.7225 -102.6567
#2007.1532153411974 -74.72270202636719 -102.65679931640625
#2007.1533237886858 -74.6662 -102.4143
#2007.1534037290714 -74.66619873046875 -102.41439819335938
#2007.1536320395737 -74.6198 -101.9188
#2007.153687182902 -74.62000274658203 -101.91840362548828

##

#1994.1985901826483 -75.061 -101.848
#2007.153687182902 -74.62000274658203 -101.91840362548828

# PIG_S mooring: -75.0546, -102.1588	
# PIG_N mooring: -74.864, -102.0989


quit()









