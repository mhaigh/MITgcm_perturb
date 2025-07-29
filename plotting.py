# plotting.py

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches

from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.basemap import Basemap

from plotting_tools import setText, getContourfLevels, maskBathyXY, maskDraftXY, makeList, doTicks, doXticks, doYticks, doTitle, doXlabels, doYlabels, doLabels, doStippling, doBox, makeBathyContour, doAntarcticInset, makeStaggeredGrid

#=============================================================================================

black = 'k'; blue = 'cornflowerblue'; red = 'indianred'; green = 'seagreen'; grey = 'gray'
cmapMean = 'plasma'; cmapAnom = 'coolwarm'; white = 'w'; cyan = 'cyan'; darkgrey = 'dimgrey'

#=============================================================================================

def quiver1byNBasemap(u, v, Xd, Yd, lat_0, lon_0, C=None, ccmap='bwr', contourf=None, contourfNlevels=9, X=None, Y=None, mesh=False, isf=None, contour=None, contourLevels=None, ls='solid', lw=0.4, cc='k', contour2=None, contourLevels2=None, ls2='solid', lw2=0.4, cc2='grey', contour3=None, contourLevels3=None, ls3='solid', lw3=0.4, cc3='w', contour4=None, contourLevels4=None, ls4='solid', lw4=0.4, cc4='w', DASH_NEG_CONTOURS=[False]*4, land=None, titles=None, fstitle=None, fontsize=14, cmap='jet', vmin=None, vmax=None, xlabel=None, ylabel=None, parallels=None, meridians=None, gridLinesWidth=0.8, text_data=None,  labelData=None, cbar=True,  cbarTicks=None, cbarLabels=None, cbarShrink=1., grid=False, scale=2, qs=0.1, qunits='m/s', labelpos='N', extend='neither', width=0.003, headwidth=5., headlength=5., headaxislength=3, qcolor='k', qlabelx=0.3, qlabely=0.03, figsize=(8,8), save=False, outpath='', outname='quiver1byN.png', show=True, dpi=200, AntarcticInsetData=None, arrows=None, plotvlines=None, maskColour='0.6'):

	#=
	#=
	
	N = len(u)
	fig = plt.figure(figsize=figsize, dpi=dpi)#, constrained_layout=True)
	gs = gridspec.GridSpec(ncols=1, nrows=N, figure=fig)

	# Basemap column

	Xd = makeList(Xd, N)
	Yd = makeList(Yd, N)

	if contourf is not None:
		X = makeList(X, N)
		Y = makeList(Y, N)

	vmin = makeList(vmin, N)
	vmax = makeList(vmax, N)
	contourfNlevels = makeList(contourfNlevels, N)
	
	AntarcticInsetData = makeList(AntarcticInsetData, N)
	arrows = makeList(arrows, N)
	plotvlines = makeList(plotvlines, N)
	
	contour = makeList(contour, N)
	contour2 = makeList(contour2, N)
	contour3 = makeList(contour3, N)
	contour4 = makeList(contour4, N)
	land = makeList(land, N)
	
	lw = makeList(lw, N)
	lw2 = makeList(lw2, N)
	lw3 = makeList(lw3, N)
	lw4 = makeList(lw4, N)	

	cmap = makeList(cmap, N)
	cbar = makeList(cbar, N)
	extend = makeList(extend, N)
	cbarTicks = makeList(cbarTicks, N)
	cbarLabels = makeList(cbarLabels, N)
	labelData = makeList(labelData, N)
		
	grid = makeList(grid, N)
	titles = makeList(titles, N)
	fstitle = makeList(fstitle, N)
			
	xlabel = makeList(xlabel, N)
	ylabel = makeList(ylabel, N)

	qcolor = makeList(qcolor, N)
	qs = makeList(qs, N)
	scale = makeList(scale, N)
	labelpos = makeList(labelpos, N)
	
	m = Basemap(llcrnrlon=X[0][0,0],llcrnrlat=Y[0][0,0],urcrnrlon=X[0][-1,-1],urcrnrlat=Y[0][-1,-1], projection='merc',lat_0=lat_0+3,lon_0=lon_0-10)
	for pi in range(N):

		ax = fig.add_subplot(gs[pi])
			
		ax.patch.set_color(maskColour)

		X0, Y0 = m(X[pi],Y[pi])
		Xd0, Yd0 = m(Xd[pi],Yd[pi])
	
		if contourf is not None:
			if mesh:
				plt.pcolormesh(X0, Y0, contourf[pi], vmin=vmin[pi], vmax=vmax[pi], cmap=cmap[pi])		
			else:
				levels = getContourfLevels(vmin[pi], vmax[pi], contourfNlevels[pi])
				plt.contourf(X0, Y0, contourf[pi], cmap=cmap[pi], levels=levels, extend=extend[pi])

			if cbar[pi]:
				cbar_ = plt.colorbar(ticks=cbarTicks[pi], shrink=cbarShrink)
				cbar_.ax.tick_params(labelsize=fontsize)
				cbar_.set_label(cbarLabels[pi], rotation=270)
				cbar_.ax.get_yaxis().labelpad = 13
				
		if contour[pi] is not None:
			if DASH_NEG_CONTOURS[0]:
				ls = np.where(np.array(contourLevels[pi]) > 0, "-", "--")
			plt.contour(X0, Y0, contour[pi], colors=cc[pi], linestyles=ls, linewidths=lw[pi], levels=contourLevels[pi],zorder=1)
		if contour2[pi] is not None:
			if DASH_NEG_CONTOURS[1]:
				ls2 = np.where(np.array(contourLevels2[pi]) > 0, "-", "--")
			plt.contour(X0, Y0, contour2[pi], colors=cc2[pi], linestyles=ls2, linewidths=lw2[pi], levels=contourLevels2[pi],zorder=1)
			
		if contour3[pi] is not None:
			if DASH_NEG_CONTOURS[2]:
				ls3 = np.where(np.array(contourLevels3[pi]) > 0, "-", "--")
			plt.contour(X0, Y0, contour3[pi], colors=cc3[pi], linestyles=ls3, linewidths=lw3[pi], levels=contourLevels3[pi], zorder=1)
		if contour4[pi] is not None:
			if DASH_NEG_CONTOURS[3]:
				ls4 = np.where(np.array(contourLevels4[pi]) > 0, "-", "--")
			plt.contour(X0, Y0, contour4[pi], colors=cc4[pi], linestyles=ls4, linewidths=lw4[pi], levels=contourLevels4[pi], zorder=14)

		if isf[pi] is not None:
			extent = [X0[0,0], X0[0,-1], -Y0[0,0], -Y0[-1,0]]
			m.imshow(1-isf[pi], cmap=plt.cm.gray, interpolation='nearest', extent=extent, zorder=13)

		if land[pi] is not None:
			extent = [X[pi][0,0], X[pi][0,-1], -Y[pi][0,0], -Y[pi][-1,0]]
			m.imshow(1-land[pi], cmap='gray', interpolation='nearest', extent=extent,);
		
		if u[pi] is not None:
			if C is not None:
				cax = ax.quiver(Xd0, Yd0, u[pi], v[pi], C[pi], cmap=ccmap, scale=scale[pi])
				plt.colorbar(cax, ax=ax)
			else:
				cax = plt.quiver(Xd0, Yd0, u[pi], v[pi], scale=scale[pi], width=width[pi], headwidth=headwidth[pi], headlength=headlength[pi], headaxislength=headaxislength[pi], color=qcolor[pi])
			ax.quiverkey(cax, qlabelx[pi], qlabely[pi], qs[pi], str(qs[pi]) + ' ' + qunits[pi], labelpos=labelpos[pi], coordinates='axes', labelcolor=qcolor[pi], fontproperties={'size':fontsize+4})
				
		doLabels(xlabel[pi], ylabel[pi], fontsize=fontsize)
		doTitle(titles[pi], fontsize=fstitle[pi])

		if text_data[pi] is not None:
			setText(ax, text_data[pi], set_invisible=False)
		
		if labelData[pi] is not None:
			for li in labelData[pi]:
				plt.scatter(li['x'][0], li['x'][1], s=1, color='r')
				plt.annotate(li['t'], li['x'])

		#ax.set_aspect('equal')

		if grid[pi]:
			plt.grid()
		
		ax.set_aspect('equal')
		if parallels[pi] is not None:
			# Latitudes
		 	m.drawparallels(parallels[pi],labels=[True,False,False,True], color='k', linewidth=1.e-6, dashes=[4,4], fontsize=fontsize)
		if meridians[pi] is not None:
			# Longitudes
			m.drawmeridians(meridians[pi],labels=[True,False,False,True], color='k', linewidth=1.e-6, dashes=[4,4], fontsize=fontsize)
			
		if plotvlines[pi] is not None:
			for vline in plotvlines[pi]:
				nvline = 40
				m.plot(vline[0]*np.ones(nvline), np.linspace(vline[1], vline[2], nvline), color=vline[3], linestyle='dashed', latlon=True, lw=vline[4])
			
		if arrows[pi] is not None:
			for arrow in arrows[pi]:
				plt.arrow(arrow[0], arrow[1], arrow[2], arrow[3], color=arrow[4], zorder=15, linewidth=3, head_width=1.5e4, head_length=1.5e4)
				
		if AntarcticInsetData[pi] is not None:
			doAntarcticInset(fig, AntarcticInsetData)

	#==

	#plt.subplots_adjust(hspace=0.05, wspace=0.05)
	#plt.subplots_adjust(left=spleft, right=spright)
	#fig.subplots_adjust(bottom=spbottom, top=sptop)

	plt.tight_layout()

	if save:
		plt.savefig(outpath + outname, bbox_inches="tight")

	if show:
		plt.show()

	plt.close()

#==

def plot1by1Basemap_Hovmoller_OLD(data, X, Y, lat_0, lon_0, h_data, h_time, h_Z,\
qx=0.3, qy=0.03, qs=0.05, qunits='m/s', scale=1, vmin=None, vmax=None, width=2.8e6, height=1.7e6, stippling=None, stipData=[0.05, 4, 4, .1], stipLim2=None, plotvlines=[], contourf=None, contourfNlevels=13, contour=None, lw=1.2, contourLevels=None, contourColours=None, contourLS=None, contourLW=None, figsize=(5,4), title='', fstitle=14, fontsize=14, mesh=True, cmap='jet', xlabel='', ylabel='', text_data=None, parallels=None, meridians=None, isf=None, land=None, outline=None, cbarTicks=None, extend='neither', AntarcticInsetData=None, maskColour='.6',\
h_cmap='coolwarm', h_mesh=False, h_xlabel=None, h_ylabel=None, h_title='', h_vmin=None, h_vmax=None, h_extend='both', h_fontsize=14, h_contourfNlevels=15, \
save=False, outpath='', outname='plot1by1Basemap_Hovmoller.png', show=True, dpi=200, height_ratios=None):

	m = Basemap(llcrnrlon=X[0,0],llcrnrlat=Y[0,0],urcrnrlon=X[-1,-1],urcrnrlat=Y[-1,-1], projection='merc',lat_0=lat_0+3,lon_0=lon_0-10)

	X,Y = m(X,Y)

	fig = plt.figure(figsize=figsize, dpi=dpi)
	
	gs = gridspec.GridSpec(ncols=1, nrows=2, figure=fig, height_ratios=height_ratios)
	ax = fig.add_subplot(gs[0])
	plt.gca().patch.set_color(maskColour)

	if mesh:
		m.pcolormesh(X, Y, data, cmap=cmap, vmin=vmin, vmax=vmax)
		cbar = m.colorbar(extend=extend, ticks=cbarTicks)
		cbar.ax.tick_params(labelsize=fontsize+1)
	else:
		levels = getContourfLevels(vmin, vmax, contourfNlevels)
		m.contourf(X, Y, data, levels=levels, cmap=cmap)
		cbar = m.colorbar(extend=extend, ticks=cbarTicks)
		cbar.ax.tick_params(labelsize=fontsize+1)
	
	if contour is not None:	
		for ci in range(len(contour)):
			cs = m.contour(X, Y, contour[ci], levels=contourLevels[ci], colors=contourColours[ci], linestyles=contourLS[ci], linewidths=contourLW[ci])
			#plt.clabel(cs, fontsize=7)
	
	if stippling is not None:
		doStippling(stippling, X=X, Y=Y, stipData=stipData, stipLim2=stipLim2)
		
	if xlabel is not None:
		plt.xlabel(xlabel, labelpad=20)
	if ylabel is not None:
		plt.ylabel(ylabel, labelpad=35)

	if title is not None:
		plt.title(title, fontsize=fstitle) 

	if text_data is not None:
		setText(ax, text_data)

	ax.set_aspect(1)
	if parallels is not None:
		# Latitudes
		m.drawparallels(parallels,labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[4,4], zorder=4, fontsize=fontsize+2)
	if meridians is not None:
		# Longitudes
		m.drawmeridians(meridians,labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[4,4], zorder=4, fontsize=fontsize+2)

	if isf is not None:
		extent = [X[0,0], X[0,-1], -Y[0,0], -Y[-1,0]]
		m.imshow(1-isf, cmap=plt.cm.gray, interpolation='nearest', extent=extent, zorder=3);
		#m.contourf(X, Y, 1-isf, cmap=plt.cm.gray, extent=extent, zorder=3, vmin=0, vmax=1);
		
	if land is not None:
		extent = [X[0,0], X[0,-1], -Y[0,0], -Y[-1,0]]
		m.imshow(1-land, cmap='gray', interpolation='nearest', extent=extent, zorder=2);
		
	if outline is not None:
		lons = outline[0]
		lats = outline[1]
		lons, lats = m(lons, lats)

		m.plot([lons[0], lons[0]], [lats[0], lats[1]], color='k', linestyle='--')
		m.plot([lons[1], lons[1]], [lats[0], lats[1]], color='k', linestyle='--')
		m.plot([lons[0], lons[1]], [lats[0], lats[0]], color='k', linestyle='--')
		m.plot([lons[0], lons[1]], [lats[1], lats[1]], color='k', linestyle='--')


	for vline in plotvlines:
		nvline = 40
		m.plot(vline[0]*np.ones(nvline), np.linspace(vline[1], vline[2], nvline), color=vline[3], linestyle='dashed', latlon=True, lw=vline[4])

	if AntarcticInsetData is not None:
		doAntarcticInset(fig, AntarcticInsetData)
		#doAntarcticInset(ax, AntarcticInsetData)
		
	#==
	
	# Now Hovmoller diagram
			
	ax = fig.add_subplot(gs[1])
	
	if h_mesh:
		plt.pcolormesh(h_time, h_Z, h_data, cmap=h_cmap)
	else:
		levels = getContourfLevels(h_vmin, h_vmax, h_contourfNlevels)
		plt.contourf(h_time, h_Z, h_data, cmap=h_cmap, levels=levels)

	if h_title is not None:
		plt.title(h_title, fontsize=fstitle)

	if h_xlabel is not None:
		plt.xlabel(h_xlabel, fontsize=h_fontsize)
	if h_ylabel is not None:
		plt.ylabel(h_ylabel, fontsize=h_fontsize)
	
	ax.tick_params(axis='both', which='major', labelsize=h_fontsize)
	ax.tick_params(axis='both', which='minor', labelsize=h_fontsize)

	plt.grid()
	plt.colorbar(extend=h_extend)
		
	#==	

	plt.tight_layout()

	if save:
		plt.savefig(outpath + outname)
		
	if show:
		plt.show()
		
#==

def plot1by1_Hovmoller(data, X, Y, h_data, h_time, h_Z,\
vmin=None, vmax=None, width=2.8e6, height=1.7e6, stippling=None, stipData=[0.05, 4, 4, .1], stipLim2=None, plotvlines=[], contourf=None, contourfNlevels=13, contour=None, lw=1.2, contourLevels=None, contourColours=None, contourLS=None, contourLW=None, figsize=(5,4), title='', fstitle=14, fontsize=14, mesh=True, cmap='jet', xlabel='', ylabel='', text_data=None, parallels=None, paraLabels=None, meridians=None, meridLabels=None, isf=None, land=None, outline=None, cbarTicks=None, extend='neither', AntarcticInsetData=None, maskColour='.6',\
h_cmap='coolwarm', h_mesh=False, h_xlabel=None, h_ylabel=None, h_title='', h_vmin=None, h_vmax=None, h_extend='both', h_fontsize=14, h_contourfNlevels=15, h_timeseries=None, h_ylim=None, h_tylim=None, h_tyticks=None,\
save=False, outpath='', outname='plot1by1Basemap_Hovmoller.png', show=True, dpi=200, height_ratios=None):


	#fig = plt.figure(figsize=figsize, dpi=dpi)
	fig,ax = plt.subplots(figsize=figsize, dpi=dpi, nrows=2, gridspec_kw=dict(height_ratios=[1.4,1]))
	
	plt.subplot(211)
	
	ax = plt.gca()
	ax.patch.set_color(maskColour)

	if mesh:
		plt.pcolormesh(X, Y, data, cmap=cmap, vmin=vmin, vmax=vmax)
		cbar = plt.colorbar(extend=extend, ticks=cbarTicks)
		cbar.ax.tick_params(labelsize=fontsize+1)
	else:
		levels = getContourfLevels(vmin, vmax, contourfNlevels)
		plt.contourf(X, Y, data, levels=levels, cmap=cmap)
		cbar = plt.colorbar(extend=extend, ticks=cbarTicks)
		cbar.ax.tick_params(labelsize=fontsize+1)
	
	if contour is not None:	
		for ci in range(len(contour)):
			cs = plt.contour(X, Y, contour[ci], levels=contourLevels[ci], colors=contourColours[ci], linestyles=contourLS[ci], linewidths=contourLW[ci])
			#plt.clabel(cs, fontsize=7)
	
	if stippling is not None:
		doStippling(stippling, X=X, Y=Y, stipData=stipData, stipLim2=stipLim2)
		
	if xlabel is not None:
		plt.xlabel(xlabel, labelpad=20)
	if ylabel is not None:
		plt.ylabel(ylabel, labelpad=35)

	plt.yticks(parallels, paraLabels, fontsize=fontsize)
	plt.xticks(meridians, meridLabels, fontsize=fontsize)
	plt.grid(linestyle='--', color='k', linewidth=0.4)
	
	if title is not None:
		plt.title(title, fontsize=fstitle) 

	if text_data is not None:
		setText(ax, text_data)

	if isf is not None:
		extent = [X[0,0], X[0,-1], Y[-1,0], Y[0,0]]
		#plt.imshow(1-np.roll(isf[:,:],0,axis=0), cmap=plt.cm.gray, interpolation='nearest', extent=extent, zorder=3, aspect='auto');
		plt.pcolormesh(X,Y,isf,cmap=plt.cm.gray,alpha=np.where(isf==1,1,0),zorder=3)
		
	if land is not None:
		extent = [X[0,0], X[0,-1], -Y[0,0], -Y[-1,0]]
		plt.imshow(1-land, cmap='gray', interpolation='nearest', extent=extent, zorder=2);
		
	if outline is not None:
		lons = outline[0]
		lats = outline[1]
		lons, lats = m(lons, lats)

		plt.plot([lons[0], lons[0]], [lats[0], lats[1]], color='k', linestyle='--')
		plt.plot([lons[1], lons[1]], [lats[0], lats[1]], color='k', linestyle='--')
		plt.plot([lons[0], lons[1]], [lats[0], lats[0]], color='k', linestyle='--')
		plt.plot([lons[0], lons[1]], [lats[1], lats[1]], color='k', linestyle='--')


	for vline in plotvlines:
		nvline = 40
		plt.plot(vline[0]*np.ones(nvline), np.linspace(vline[1], vline[2], nvline), color=vline[3], linestyle='dashed', latlon=True, lw=vline[4])

	if AntarcticInsetData is not None:
		doAntarcticInset(fig, AntarcticInsetData)
		#doAntarcticInset(ax, AntarcticInsetData)
	
	ax.tick_params(axis='both', which='major', labelsize=h_fontsize)
	ax.tick_params(axis='both', which='minor', labelsize=h_fontsize)
	
	plt.grid()
	
	#==
	
	# Now Hovmoller diagram
			
	plt.subplot(212)
	ax = plt.gca()
	
	if h_mesh:
		plt.pcolormesh(h_time, h_Z, h_data, cmap=h_cmap)
	else:
		levels = getContourfLevels(h_vmin, h_vmax, h_contourfNlevels)
		plt.contourf(h_time, h_Z, h_data, cmap=h_cmap, levels=levels, extend=h_extend)
	cbar = plt.colorbar()
	cbar.ax.tick_params(labelsize=fontsize+1)
		
	if h_title is not None:
		plt.title(h_title, fontsize=fstitle)

	if h_xlabel is not None:
		plt.xlabel(h_xlabel, fontsize=h_fontsize)
	if h_ylabel is not None:
		plt.ylabel(h_ylabel, fontsize=h_fontsize)
	
	ax.tick_params(axis='both', which='major', labelsize=h_fontsize)
	ax.tick_params(axis='both', which='minor', labelsize=h_fontsize)
	plt.ylim(h_ylim)
	plt.grid()

	if h_timeseries is not None:
		ax2 = ax.twinx()
		ax2.plot(h_time, h_timeseries, color=green)
	plt.yticks(h_tyticks, fontsize=fontsize, color=green)
	plt.ylim(h_tylim)

		
	#==	

	plt.tight_layout()

	if save:
		plt.savefig(outpath + outname)
		
	if show:
		plt.show()

#==

def plot1by1_timeseries(data, X, Y,\
vmin=None, vmax=None, width=2.8e6, height=1.7e6, stippling=None, stipData=[0.05, 4, 4, .1], stipLim2=None, plotvlines=[], contourf=None, contourfNlevels=13, contour=None, lw=1.2, contourLevels=None, contourColours=None, contourLS=None, contourLW=None, figsize=(5,4), title='', fstitle=14, fontsize=14, mesh=True, cmap='jet', xlabel='', ylabel='', text_data=None, parallels=None, paraLabels=None, meridians=None, meridLabels=None, isf=None, land=None, outline=None, cbarTicks=None, extend='neither', AntarcticInsetData=None, maskColour='.6',\
ts_time=None, ts_data=None, ts_labels=None, ts_lines=None, ts_title=None, ts_xlabel=None, ts_ylabel=None, ts_xlim=None, ts_ylim=None, ts_fslegend=9,\
save=False, outpath='', outname='plot1by1Basemap_timeseries.png', show=True, dpi=200, height_ratios=None):


	#fig = plt.figure(figsize=figsize, dpi=dpi)
	fig,ax = plt.subplots(figsize=figsize, dpi=dpi, nrows=2, gridspec_kw=dict(height_ratios=[1.4,1]))
	
	plt.subplot(211)
	
	ax = plt.gca()
	ax.patch.set_color(maskColour)

	if mesh:
		plt.pcolormesh(X, Y, data, cmap=cmap, vmin=vmin, vmax=vmax)
		cbar = plt.colorbar(extend=extend, ticks=cbarTicks)
		cbar.ax.tick_params(labelsize=fontsize+1)
	else:
		levels = getContourfLevels(vmin, vmax, contourfNlevels)
		plt.contourf(X, Y, data, levels=levels, cmap=cmap)
		cbar = plt.colorbar(extend=extend, ticks=cbarTicks)
		cbar.ax.tick_params(labelsize=fontsize+1)
	
	if contour is not None:	
		for ci in range(len(contour)):
			cs = plt.contour(X, Y, contour[ci], levels=contourLevels[ci], colors=contourColours[ci], linestyles=contourLS[ci], linewidths=contourLW[ci])
			#plt.clabel(cs, fontsize=7)
	
	if stippling is not None:
		doStippling(stippling, X=X, Y=Y, stipData=stipData, stipLim2=stipLim2)
		
	if xlabel is not None:
		plt.xlabel(xlabel, labelpad=20)
	if ylabel is not None:
		plt.ylabel(ylabel, labelpad=35)

	plt.yticks(parallels, paraLabels, fontsize=fontsize)
	plt.xticks(meridians, meridLabels, fontsize=fontsize)
	plt.grid(linestyle='--', color='k', linewidth=0.4)
	
	if title is not None:
		plt.title(title, fontsize=fstitle) 

	if text_data is not None:
		setText(ax, text_data)

	if isf is not None:
		extent = [X[0,0], X[0,-1], Y[-1,0], Y[0,0]]
		#plt.imshow(1-np.roll(isf[:,:],0,axis=0), cmap=plt.cm.gray, interpolation='nearest', extent=extent, zorder=3, aspect='auto');
		plt.pcolormesh(X,Y,isf,cmap=plt.cm.gray,alpha=np.where(isf==1,1,0),zorder=3)
		
	if land is not None:
		extent = [X[0,0], X[0,-1], -Y[0,0], -Y[-1,0]]
		plt.imshow(1-land, cmap='gray', interpolation='nearest', extent=extent, zorder=2);
		
	if outline is not None:
		lons = outline[0]
		lats = outline[1]
		lons, lats = m(lons, lats)

		plt.plot([lons[0], lons[0]], [lats[0], lats[1]], color='k', linestyle='--')
		plt.plot([lons[1], lons[1]], [lats[0], lats[1]], color='k', linestyle='--')
		plt.plot([lons[0], lons[1]], [lats[0], lats[0]], color='k', linestyle='--')
		plt.plot([lons[0], lons[1]], [lats[1], lats[1]], color='k', linestyle='--')


	for vline in plotvlines:
		nvline = 40
		plt.plot(vline[0]*np.ones(nvline), np.linspace(vline[1], vline[2], nvline), color=vline[3], linestyle='dashed', latlon=True, lw=vline[4])

	if AntarcticInsetData is not None:
		doAntarcticInset(fig, AntarcticInsetData)
		#doAntarcticInset(ax, AntarcticInsetData)
	
	ax.tick_params(axis='both', which='major', labelsize=fontsize)
	ax.tick_params(axis='both', which='minor', labelsize=fontsize)
	
	plt.grid()
	
	#==
	
	# Now timeseries
	
	plt.subplot(212)
	ax = plt.gca()

	for ti, timeseries in enumerate(ts_data):
		plt.plot(ts_time[ti], timeseries, color=ts_lines[ti][0], linestyle=ts_lines[ti][1], label=ts_labels[ti]) 
		
	plt.xlim(ts_xlim)
	plt.ylim(ts_ylim)
	plt.title(ts_title, fontsize=fstitle)
	#ax.tick_params(axis='both', which='minor', labelsize=fontsize)
	plt.xticks(fontsize=fontsize)
	plt.yticks(fontsize=fontsize)
	plt.xlabel(ts_xlabel, fontsize=fontsize)
	plt.ylabel(ts_ylabel, fontsize=fontsize)
	plt.grid(color='k', linestyle='dashed', linewidth=0.5)
	plt.legend(prop={'size': ts_fslegend}, ncol=1)

	
	#==	

	plt.tight_layout()

	if save:
		plt.savefig(outpath + outname)
		
	if show:
		plt.show()
		
		
#==

def quiverMbyN_Basemap(u, v, Xd, Yd, lat_0, lon_0, C=None, ccmap='bwr', contourf=None, cbarTicks=None, contourfNlevels=9, X=None, Y=None, mesh=False, isf=None, contour=None, contourLevels=None, ls='solid', lw=0.4, cc='k', contour2=None, contourLevels2=None, ls2='solid', lw2=0.4, cc2='grey', contour3=None, contourLevels3=None, ls3='solid', lw3=0.4, cc3='w', contour4=None, contourLevels4=None, ls4='solid', lw4=0.4, cc4='w', DASH_NEG_CONTOURS=[False]*4, figsize=None, titles=None, fstitle=None, fontsize=14, cmap='jet', vmin=None, vmax=None, xlabel=None, ylabel=None, parallels=None, meridians=None, save=False, outpath='', outname='quiver1byN.png', show=True, dpi=200, text_data=None, width_ratios=1, labelData=None, cbar=True, cbarShrink=1., grid=True, scale=2, qs=0.1, qunits='m/s', labelpos='E', landscape=True, extend='neither', width=0.003, headwidth=5., headlength=5., headaxislength=3, qcolor='k', qlabelx=0.3, qlabely=0.03, w_pad=0, h_pad=0, drawCoastlines=False, landIce=None, XlandIce=None, YlandIce=None,\
t_=None):

	N = len(u)
	M = len(u[0])
	
	Xd = makeList(Xd, M, N)
	Yd = makeList(Yd, M, N)

	if contourf is not None:
		X = makeList(X, M, N)
		Y = makeList(Y, M, N)

	vmin = makeList(vmin, M, N)
	vmax = makeList(vmax, M, N)

	contour = makeList(contour, M, N)
	contour2 = makeList(contour2, M, N)
	contour3 = makeList(contour3, M, N)
	contour4 = makeList(contour4, M, N)
	
	isf = makeList(isf, M, N)
	parallels = makeList(parallels, M, N)
	meridians = makeList(meridians, M, N)
	drawCoastlines = makeList(drawCoastlines, M, N)
	landIce = makeList(landIce, M, N)
	
	cmap = makeList(cmap, M, N)
	cbar = makeList(cbar, M, N)
	extend = makeList(extend, M, N)
	cbarTicks = makeList(cbarTicks, M, N)
	contourfNlevels = makeList(contourfNlevels, M, N)
	
	grid = makeList(grid, M, N)
	titles = makeList(titles, M, N)
	if fstitle is None:
		fstitle = makeList(fontsize, M, N)
		
	width_ratios = makeList(width_ratios, M)

	xlabel = makeList(xlabel, M, N)
	ylabel = makeList(ylabel, M, N)

	qcolor = makeList(qcolor, M, N)
	qs = makeList(qs, M, N)
	scale = makeList(scale, M, N)
	labelpos = makeList(labelpos, M, N)
	
	m = Basemap(llcrnrlon=X[0][0][0,0],llcrnrlat=Y[0][0][0,0],urcrnrlon=X[0][0][-1,-1],urcrnrlat=Y[0][0][-1,-1], projection='merc',lat_0=lat_0+3,lon_0=lon_0-10)
	
	fig = plt.figure(figsize=figsize, dpi=dpi)#, constrained_layout=True)
	gs = gridspec.GridSpec(ncols=M, nrows=N, figure=fig, width_ratios=width_ratios)
	
	XlandIce, YlandIce = m(XlandIce,YlandIce)
				
	for col in range(M):
		for row in range(N):

			ax = fig.add_subplot(gs[row,col])
			ax.patch.set_color('.6')

			X0, Y0 = m(X[row][col],Y[row][col])
			Xd0, Yd0 = m(Xd[row][col],Yd[row][col])

			if contourf is not None:
				if mesh:
					m.pcolormesh(X0, Y0, contourf[row][col], vmin=vmin[row][col], vmax=vmax[row][col], cmap=cmap[row][col])		
				else:
					levels = getContourfLevels(vmin[row][col], vmax[row][col], contourfNlevels[row][col])
					m.contourf(X0, Y0, contourf[row][col], cmap=cmap[row][col], levels=levels, extend=extend[row][col])

				if cbar[row][col]:
					cbar_ = plt.colorbar(ticks=cbarTicks[row][col], shrink=cbarShrink)
					cbar_.ax.tick_params(labelsize=fontsize)
			
			if contour[row][col] is not None:
				if DASH_NEG_CONTOURS[0]:
					ls = np.where(np.array(contourLevels[row][col]) > 0, "-", "--")
				m.contour(XlandIce, YlandIce, contour[row][col], colors=cc, linestyles=ls, linewidths=lw, levels=contourLevels[row][col])
			if contour2[row][col] is not None:
				if DASH_NEG_CONTOURS[1]:
					ls2 = np.where(np.array(contourLevels2[row][col]) > 0, "-", "--")
				m.contour(XlandIce, YlandIce, contour2[row][col], colors=cc2, linestyles=ls2, linewidths=lw2, levels=contourLevels2[row][col],zorder=1)
			if contour3[row][col] is not None:
				if DASH_NEG_CONTOURS[2]:
					ls3 = np.where(np.array(contourLevels3[row][col]) > 0, "-", "--")
				m.contour(X0, Y0, contour3[row][col], colors=cc3, linestyles=ls3, linewidths=lw3, levels=contourLevels3[row][col], zorder=1)
			if contour4[row][col] is not None:
				if DASH_NEG_CONTOURS[3]:
					ls4 = np.where(np.array(contourLevels4[row][col]) > 0, "-", "--")
				m.contour(X0, Y0, contour4[row][col], colors=cc4, linestyles=ls4, linewidths=lw4, levels=contourLevels4[row][col], zorder=1)

			if isf[row][col] is not None:
				extent = [X0[0,0], X0[0,-1], -Y0[0,0], -Y0[-1,0]]
				m.imshow(1-isf[row][col], cmap=plt.cm.gray, interpolation='nearest', extent=extent, zorder=13)
			
			if landIce[row][col] is not None:
				m.pcolormesh(XlandIce, YlandIce, landIce[row][col], alpha=np.where(landIce[row][col]<1,1,0), cmap=plt.cm.gray,zorder=2) 
			#ax.set_aspect('equal')
			
			if u[row][col] is not None:
				if C is not None:
					cax = ax.quiver(Xd0, Yd0, u[row][col], v[row][col], C[row][col], cmap=ccmap, scale=scale[row][col])
					if width_ratios is None:
						plt.colorbar(cax, ax=ax)
				else:
					cax = plt.quiver(Xd0, Yd0, u[row][col], v[row][col], scale=scale[row][col], width=width, headwidth=headwidth, headlength=headlength, headaxislength=headaxislength, color=qcolor[row][col], zorder=4)
				ax.quiverkey(cax, qlabelx, qlabely, qs[row][col], str(qs[row][col]) + ' ' + qunits, labelpos=labelpos[row][col], fontproperties={'size':fontsize+2}, coordinates='axes', labelcolor=qcolor[row][col])
					
			doLabels(xlabel[row][col], ylabel[row][col], fontsize=fontsize)
			doTitle(titles[row][col], fontsize=fstitle[row][col])

			if text_data is not None:
				setText(ax, text_data[row][col], set_invisible=False)
			
			if labelData is not None:
				for li in labelData[row][col]:
					plt.scatter(li['x'][0], li['x'][1], s=1, color='r')
					plt.annotate(li['t'], li['x'])


			if grid[row][col]:
				plt.grid()
			
			ax.set_aspect('equal')
			if parallels[row][col] is not None:
				# Latitudes
				m.drawparallels(parallels[row][col],labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[4,4], fontsize=fontsize, zorder=3)
			if meridians[row][col] is not None:
				# Longitudes
				m.drawmeridians(meridians[row][col],labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[4,4], fontsize=fontsize, zorder=3)

			if drawCoastlines[row][col]:
				m.fillcontinents()
				m.drawcoastlines()

		#==
	
	#fig.subplots_adjust(wspace=-1, hspace=0)

	plt.tight_layout(w_pad=w_pad, h_pad=h_pad)

	if save:
		plt.savefig(outpath + outname, bbox_inches="tight")

	if show:
		plt.show()

	plt.close()

#==

def quiverMbyN_Basemap_timeseries(u, v, Xd, Yd, lat_0, lon_0, C=None, ccmap='bwr', contourf=None, cbarTicks=None, contourfNlevels=9, X=None, Y=None, mesh=False, isf=None, contour=None, contourLevels=None, ls='solid', lw=0.4, cc='k', contour2=None, contourLevels2=None, ls2='solid', lw2=0.4, cc2='grey', contour3=None, contourLevels3=None, ls3='solid', lw3=0.4, cc3='w', contour4=None, contourLevels4=None, ls4='solid', lw4=0.4, cc4='w', DASH_NEG_CONTOURS=[False]*4, figsize=None, titles=None, fstitle=None, fontsize=14, cmap='jet', vmin=None, vmax=None, xlabel=None, ylabel=None, parallels=None, meridians=None, save=False, outpath='', outname='quiver1byN.png', show=True, dpi=200, text_data=None, height_ratios=1, width_ratios=1, labelData=None, cbar=True, cbarShrink=1., grid=True, scale=2, qs=0.1, qunits='m/s', labelpos='E', landscape=True, extend='neither', width=0.003, headwidth=5., headlength=5., headaxislength=3, qcolor='k', qlabelx=0.3, qlabely=0.03, w_pad=0, h_pad=0, drawCoastlines=False, landIce=None, XlandIce=None, YlandIce=None,\
ts_time=None, ts_data=None, ts_labels=None, ts_lines=None, ts_title=None, ts_xlabel=None, ts_ylabel=None, ts_xlim=None, ts_ylim=None, ts_fslegend=9, ts_arrows=[]):

	N = len(u)
	M = len(u[0])
	
	Xd = makeList(Xd, M, N)
	Yd = makeList(Yd, M, N)

	if contourf is not None:
		X = makeList(X, M, N)
		Y = makeList(Y, M, N)

	vmin = makeList(vmin, M, N)
	vmax = makeList(vmax, M, N)

	contour = makeList(contour, M, N)
	contour2 = makeList(contour2, M, N)
	contour3 = makeList(contour3, M, N)
	contour4 = makeList(contour4, M, N)
	
	isf = makeList(isf, M, N)
	parallels = makeList(parallels, M, N)
	meridians = makeList(meridians, M, N)
	drawCoastlines = makeList(drawCoastlines, M, N)
	landIce = makeList(landIce, M, N)
	
	cmap = makeList(cmap, M, N)
	cbar = makeList(cbar, M, N)
	extend = makeList(extend, M, N)
	cbarTicks = makeList(cbarTicks, M, N)
	contourfNlevels = makeList(contourfNlevels, M, N)
	
	grid = makeList(grid, M, N)
	titles = makeList(titles, M, N)
	if fstitle is None:
		fstitle = makeList(fontsize, M, N)
		
	width_ratios = makeList(width_ratios, M)

	xlabel = makeList(xlabel, M, N)
	ylabel = makeList(ylabel, M, N)

	qcolor = makeList(qcolor, M, N)
	qs = makeList(qs, M, N)
	scale = makeList(scale, M, N)
	labelpos = makeList(labelpos, M, N)
	
	m = Basemap(llcrnrlon=X[0][0][0,0],llcrnrlat=Y[0][0][0,0],urcrnrlon=X[0][0][-1,-1],urcrnrlat=Y[0][0][-1,-1], projection='merc',lat_0=lat_0+3,lon_0=lon_0-10)
	
	fig = plt.figure(figsize=figsize, dpi=dpi)#, constrained_layout=True)
	#gs = gridspec.GridSpec(ncols=M, nrows=N, figure=fig, width_ratios=width_ratios)
	gs = gridspec.GridSpec(ncols=1, nrows=N+1, figure=fig, height_ratios=height_ratios)
	
	XlandIce, YlandIce = m(XlandIce,YlandIce)
				
	for col in range(M):
		for row in range(N):

			ax = fig.add_subplot(gs[row,col])
			ax.patch.set_color('.6')

			X0, Y0 = m(X[row][col],Y[row][col])
			Xd0, Yd0 = m(Xd[row][col],Yd[row][col])

			if contourf is not None:
				if mesh:
					m.pcolormesh(X0, Y0, contourf[row][col], vmin=vmin[row][col], vmax=vmax[row][col], cmap=cmap[row][col])		
				else:
					levels = getContourfLevels(vmin[row][col], vmax[row][col], contourfNlevels[row][col])
					m.contourf(X0, Y0, contourf[row][col], cmap=cmap[row][col], levels=levels, extend=extend[row][col])

				if cbar[row][col]:
					cbar_ = plt.colorbar(ticks=cbarTicks[row][col], shrink=cbarShrink)
					cbar_.ax.tick_params(labelsize=fontsize)
			
			if contour[row][col] is not None:
				if DASH_NEG_CONTOURS[0]:
					ls = np.where(np.array(contourLevels[row][col]) > 0, "-", "--")
				m.contour(XlandIce, YlandIce, contour[row][col], colors=cc, linestyles=ls, linewidths=lw, levels=contourLevels[row][col])
			if contour2[row][col] is not None:
				if DASH_NEG_CONTOURS[1]:
					ls2 = np.where(np.array(contourLevels2[row][col]) > 0, "-", "--")
				m.contour(XlandIce, YlandIce, contour2[row][col], colors=cc2, linestyles=ls2, linewidths=lw2, levels=contourLevels2[row][col],zorder=1)
			if contour3[row][col] is not None:
				if DASH_NEG_CONTOURS[2]:
					ls3 = np.where(np.array(contourLevels3[row][col]) > 0, "-", "--")
				m.contour(X0, Y0, contour3[row][col], colors=cc3, linestyles=ls3, linewidths=lw3, levels=contourLevels3[row][col], zorder=1)
			if contour4[row][col] is not None:
				if DASH_NEG_CONTOURS[3]:
					ls4 = np.where(np.array(contourLevels4[row][col]) > 0, "-", "--")
				m.contour(X0, Y0, contour4[row][col], colors=cc4, linestyles=ls4, linewidths=lw4, levels=contourLevels4[row][col], zorder=1)

			if isf[row][col] is not None:
				extent = [X0[0,0], X0[0,-1], -Y0[0,0], -Y0[-1,0]]
				m.imshow(1-isf[row][col], cmap=plt.cm.gray, interpolation='nearest', extent=extent, zorder=13)
			
			if landIce[row][col] is not None:
				m.pcolormesh(XlandIce, YlandIce, landIce[row][col], alpha=np.where(landIce[row][col]<1,1,0), cmap=plt.cm.gray,zorder=2) 
			#ax.set_aspect('equal')
			
			if u[row][col] is not None:
				if C is not None:
					cax = ax.quiver(Xd0, Yd0, u[row][col], v[row][col], C[row][col], cmap=ccmap, scale=scale[row][col])
					if width_ratios is None:
						plt.colorbar(cax, ax=ax)
				else:
					cax = plt.quiver(Xd0, Yd0, u[row][col], v[row][col], scale=scale[row][col], width=width, headwidth=headwidth, headlength=headlength, headaxislength=headaxislength, color=qcolor[row][col], zorder=4)
				ax.quiverkey(cax, qlabelx, qlabely, qs[row][col], str(qs[row][col]) + ' ' + qunits, labelpos=labelpos[row][col], fontproperties={'size':fontsize+2}, coordinates='axes', labelcolor=qcolor[row][col])
					
			doLabels(xlabel[row][col], ylabel[row][col], fontsize=fontsize)
			doTitle(titles[row][col], fontsize=fstitle[row][col])

			if text_data is not None:
				setText(ax, text_data[row][col], set_invisible=False)
			
			if labelData is not None:
				for li in labelData[row][col]:
					plt.scatter(li['x'][0], li['x'][1], s=1, color='r')
					plt.annotate(li['t'], li['x'])


			if grid[row][col]:
				plt.grid()
			
			ax.set_aspect('equal')
			if parallels[row][col] is not None:
				# Latitudes
				m.drawparallels(parallels[row][col],labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[4,4], fontsize=fontsize, zorder=3)
			if meridians[row][col] is not None:
				# Longitudes
				m.drawmeridians(meridians[row][col],labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[4,4], fontsize=fontsize, zorder=3)

			if drawCoastlines[row][col]:
				m.fillcontinents()
				m.drawcoastlines()

		#==
	
	#==
	
	# Now timeseries.
	ax = fig.add_subplot(gs[2])

	for ti, timeseries in enumerate(ts_data):
		plt.plot(ts_time[ti], timeseries, color=ts_lines[ti][0], linestyle=ts_lines[ti][1], label=ts_labels[ti]) 
		
	plt.xlim(ts_xlim)
	plt.ylim(ts_ylim)
	plt.title(ts_title, fontsize=fstitle[0][0])
	#ax.tick_params(axis='both', which='minor', labelsize=fontsize)
	plt.xticks(fontsize=fontsize)
	plt.yticks(fontsize=fontsize)
	plt.xlabel(ts_xlabel, fontsize=fontsize)
	plt.ylabel(ts_ylabel, fontsize=fontsize)
	plt.grid(color='k', linestyle='dashed', linewidth=0.5)
	plt.legend(prop={'size': ts_fslegend}, ncol=1)

	for arrow in ts_arrows:
		if arrow[2] == 0:
			plt.arrow(arrow[0], arrow[1], arrow[2], arrow[3], color=arrow[4], zorder=15, linewidth=3, head_width=.5, head_length=6)	
		else:
			plt.arrow(arrow[0], arrow[1], arrow[2], arrow[3], color=arrow[4], zorder=15, linewidth=3, head_width=6., head_length=.5)	
			
	#==
	#fig.subplots_adjust(wspace=-1, hspace=0)

	plt.tight_layout(w_pad=w_pad, h_pad=h_pad)

	if save:
		plt.savefig(outpath + outname, bbox_inches="tight")

	if show:
		plt.show()

	plt.close()


#==

def plot1by1(data, X=None, Y=None, contour=None, contourLevels=None, ls='solid', lw=0.4, cc='k', contour2=None, contourLevels2=None, ls2='solid', lw2=0.4, cc2='grey', contour3=None, contourLevels3=None, ls3='solid', lw3=0.4, cc3='w', figsize=(5,4), title=None, fontsize=14, mesh=False, cmap='jet', vmin=None, vmax=None, text_data=None, xlabel=None, ylabel=None, grid=True, contourfNlevels=9, save=True, outpath='', outname='plot1by1.png', show=False, dpi=200, hlines=None, xmin=[0], xmax=[1], vlines=None, xlim=None, ylim=None, stippling=None, stipData=[0.05, 4, 4, .1], box=None, patchColour='.25', isf=None, isfCmap=plt.cm.gray):
	
	
	if contourLevels2 == None:
		contourLevels2 = contourLevels
	if contourLevels3 == None:
		contourLevels3 = contourLevels
		
	fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

	plt.subplot(111)
	plt.gca().patch.set_color(patchColour)

	if mesh:
		if X is not None and Y is not None:
			cax = plt.pcolormesh(X, Y, data, cmap=cmap, vmin=vmin, vmax=vmax)
		else: 
			cax = plt.pcolormesh(data, cmap=cmap, vmin=vmin, vmax=vmax)
	else:
		levels = getContourfLevels(vmin, vmax, contourfNlevels)
		if X is not None and Y is not None:
			cax = plt.contourf(X, Y, data, cmap=cmap, levels=levels)
		else: 
			cax = plt.contourf(data, cmap=cmap, levels=levels)

	if contour is not None:
		plt.contour(X, Y, contour, colors=cc, linestyles=ls, linewidths=lw, levels=contourLevels)
	if contour2 is not None:
		plt.contour(X, Y, contour2, colors=cc2, linestyles=ls2, linewidths=lw2, levels=contourLevels2)
	if contour3 is not None:
		plt.contour(X, Y, contour3, colors=cc3, linestyles=ls3, linewidths=lw3, levels=contourLevels3)			
	
	if xlabel is not None:
		plt.xlabel(xlabel, fontsize=fontsize)
	if ylabel is not None:
		plt.ylabel(ylabel, fontsize=fontsize)
	
	if grid:
		plt.grid()

	if text_data is not None:
		setText(plt.gca(), text_data, set_invisible=False)

	plt.colorbar(cax, ax=ax)
	
	if title is not None:
		plt.title(title, fontsize=fontsize)

	#if yline is not None:
	#	plt.plot(yline, Y, color='k', linewidth=1.2)
	#	plt.axvline(x=X[X.shape[0]//2-1], color='k', linewidth=1.0, linestyle='--')

	if vlines is not None:
		for i in range(len(vlines)):
			plt.axvline(x=vlines[i], color='k', linewidth=1.0, linestyle='--')
	
	if hlines is not None:
		for i in range(len(hlines)):
			plt.axhline(y=hlines[i], xmin=xmin[i], xmax=xmax[i], color='k', linewidth=1.0, linestyle='--')
			
	if xlim is not None:
		plt.xlim(xlim)
	if ylim is not None:
		plt.ylim(ylim)
		
	if stippling is not None:
		doStippling(stippling, X=X, Y=Y, stipData=stipData)

	if isf is not None:
		extent = [X[0,0], X[0,-1], Y[0,0], Y[-1,0]]
		plt.imshow(1-isf, cmap=isfCmap, interpolation='nearest', extent=extent, zorder=13)
			
	doBox(box, X, Y)
		
	#==
	
	plt.tight_layout()
	if save:
		plt.savefig(outpath + outname)
		
	if show:
		plt.show()

#==

def plotMbyN(data, X=None, Y=None, figsize=(9,8), \
	titles=None, suptitle=None, insetTitle=None, \
	insetTitleX=0.05, insetTitleY=0.95, fontsize=12, insetTitleFontsize=10, \
	isf=None, mesh=False, cmap='coolwarm', vmin=None, vmax=None, cbar=False, \
	xlim=None, ylim=None, cbarShared=True, cbarSharedData=[[0.15, 0.05, 0.7, 0.015], None], \
	text_data=None, xlabel=None, ylabel=None, contourfNlevels=13, \
	save=False, outpath='', outname='plotMbyN.png', show=True, dpi=200, \
	width_ratios=None, xticks=None, yticks=None, xticksvis=False, yticksvis=False, xticklabels=None, yticklabels=None, \
	contour=None, contourLevels=None, contour2=None, contour2Levels=None, \
	DASH_NEG_CONTOURS=False, hlines=None, lines=None, linesYlim=None, \
	hspace=0.05, wspace=0.05, spleft=0.05, spright=0.95, sptop=0.95, spbottom=0.1, \
	# Quiver params
	u=None, v=None, Xd=None, Yd=None, \
	scale=2, qs=0.1, qunits='m/s', labelpos='E',  width=0.003, headwidth=5., \
	headlength=5., headaxislength=3, qcolor='k', qlabelx=0.3, qlabely=0.03, \
	# Timeseries params
	ts_time=None, ts_data=None, ts_labels=None, ts_lines=None, ts_title=None, \
	ts_xlabel=None, ts_ylabel=None, ts_xlim=None, ts_ylim=None, \
	ts_xticks=None, ts_yticks=None, ts_insetTitle=None, \
	ts_legendLoc=None, ts_legendFontsize=8, ts_legendNcol=1):

	if not isinstance(data, list):
		plot1by1(data, X=X, Y=Y, mesh=mesh, vmin=vmin, vmax=vmax)
		return
		
	N = len(data)
	M = len(data[0])
	
	if X is None:
		nX = data[0][0].shape[1]
		X = np.linspace(1, nX, nX)
	if Y is None:
		nY = data[0][0].shape[0]
		Y = np.linspace(1, nY, nY)	
	
	X = makeList(X, M, N)
	Y = makeList(Y, M, N)

	Xd = makeList(Xd, M, N)
	Yd = makeList(Yd, M, N)
	qs = makeList(qs, M, N)
	scale = makeList(scale, M, N)
	labelpos = makeList(labelpos, M, N)	

	hlines = makeList(hlines, M, N)
	
	vmin = makeList(vmin, M, N)
	vmax = makeList(vmax, M, N)	
	xlim = makeList(xlim, M, N, FORCE_MAKE_LIST=isinstance(xlim, tuple))
	ylim = makeList(ylim, M, N, FORCE_MAKE_LIST=isinstance(xlim, tuple))
		
	cbar = makeList(cbar, M, N)
	cmap = makeList(cmap, M, N)
	text_data = makeList(text_data, M, N)
	
	titles = makeList(titles, M, N)
	insetTitle = makeList(insetTitle, M, N)
	xlabel = makeList(xlabel, M, N)
	ylabel = makeList(ylabel, M, N)
	
	lines = makeList(lines, M, N)
	linesYlim = makeList(linesYlim, M, N)
	
	contour = makeList(contour, M, N)
	contourLevels = makeList(contourLevels, M, N, FORCE_MAKE_LIST=True)
	if contour2 is None:
		contour2 = makeList(contour2, M, N)		
	
	if u is None:
		u = makeList(u, M, N)
	if v is None:
		v = makeList(v, M, N)
	
	if width_ratios is None:
		width_ratios = [1]*M

	#==

	fig = plt.figure(figsize=figsize, dpi=dpi)
	#gs = gridspec.GridSpec(ncols=M, nrows=N+1, figure=fig, width_ratios=width_ratios, height_ratios=[1,1,1,1,0.05])
	gs = gridspec.GridSpec(ncols=M, nrows=N, figure=fig, width_ratios=width_ratios, height_ratios=[1]*N)
	
	for col in range(M):
		for row in range(N):

			ax = fig.add_subplot(gs[row,col])
			ax.patch.set_color('.6')

			if mesh:
				im = plt.pcolormesh(X[row][col], Y[row][col], data[row][col], vmin=vmin[row][col], vmax=vmax[row][col], cmap=cmap[row][col])		
			else:
				levels = getContourfLevels(vmin[row][col], vmax[row][col], contourfNlevels)
				im = plt.contourf(X[row][col], Y[row][col], data[row][col], cmap=cmap[row][col], levels=levels, extend='both')
			if cbar[row][col]:
				plt.colorbar()
				
			if contour[row][col] is not None:
				if DASH_NEG_CONTOURS:
					ls = np.where(np.array(contourLevels[row][col])>0,"-","--")
				else:
					ls = 'solid' 
				plt.contour(X[row][col], Y[row][col], contour[row][col], linewidths=.8, colors='k', levels=contourLevels[row][col], linestyles=ls)

			if contour2[row][col] is not None:
				plt.contour(X[row][col], Y[row][col], contour2[row][col], linewidths=.8, colors='k', levels=contour2Levels[row][col], linestyles='dashed')				

			# Quiver
			if u[row][col] is not None:
				cax = plt.quiver(Xd[row][col], Yd[row][col], u[row][col], v[row][col], scale=scale[row][col], width=width, headwidth=headwidth, headlength=headlength, headaxislength=headaxislength, color='k', zorder=4)
				ax.quiverkey(cax, qlabelx, qlabely, qs[row][col], str(qs[row][col]) + ' ' + qunits, labelpos=labelpos[row][col], fontproperties={'size':fontsize+2}, coordinates='axes', labelcolor='k')
			# End quiver
			
			if row == N-1:
				doXlabels(xlabel[row][col], fontsize=fontsize)
				doXticks(xticks, True, fontsize=10)
				if xticklabels is not None:
					ax.set_xticklabels(xticklabels)
			else:
				doXticks([], False)
			if col == 0:
				doYlabels(ylabel[row][col], fontsize=fontsize)
				doYticks(yticks, True)
				if yticklabels is not None:
					ax.set_yticklabels(yticklabels)
			else:
				doYticks([], False, fontsize=10)
			
			doTitle(titles[row][col], fontsize=fontsize)
			
			if xlim[row][col] is not None:
				plt.xlim(xlim[row][col][0], xlim[row][col][1]) 
			if ylim[row][col] is not None:
				plt.ylim(ylim[row][col][0], ylim[row][col][1]) 
				
			if hlines[row][col] is not None:
				plt.axhline(y=hlines[row][col], color='k', linewidth=0.7, linestyle='--')
		
			if text_data[row][col] is not None:
				setText(ax, text_data[row][col])
				
			if lines[row][col] is not None:
				ax2 = ax.twinx()
				ax2.plot(X[row][col], lines[row][col], color='k')
				if linesYlim[row][col] is not None:
					plt.ylim(linesYlim[row][col]) 

			if isf is not None:
				#extent = [X[row][col][0,0], X[row][col][0,-1], -Y[row][col][0,0], -Y[row][col][-1,0]]
				#plt.imshow(1-isf, cmap=plt.cm.gray, interpolation='nearest', extent=extent, zorder=1)
				plt.pcolormesh(X[row][col], Y[row][col], isf, cmap=plt.cm.gray, alpha=1)

			# Inset title
			if insetTitle[row][col] is not None:
				props = dict(boxstyle='round', facecolor='white', alpha=1.0)
				ax.text(insetTitleX, insetTitleY, insetTitle[row][col], zorder=10, transform=ax.transAxes, fontsize=insetTitleFontsize, verticalalignment='top', bbox=props)
					
	#==

	if suptitle is not None:
		plt.suptitle(suptitle, fontsize=fontsize+2)

	if cbarShared:
		if ts_data is not None:
			cbarSharedData = [[0.15, 0.233, .7, 0.015], None]
		plt.subplots_adjust(hspace=hspace, wspace=wspace)
		plt.subplots_adjust(left=spleft, right=spright)
		fig.subplots_adjust(bottom=spbottom, top=sptop)
		cbar_ax = fig.add_axes(cbarSharedData[0])
		cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
		if cbarSharedData[1] is not None:
			cbar.set_label(cbarSharedData[1], rotation=270, labelpad=20, fontsize=fontsize)
		

	#==

	if ts_data is not None:
		ax = fig.add_axes([spleft,0.025,1.-2*spleft,0.18])
		for ti, timeseries in enumerate(ts_data):
			plt.plot(ts_time[ti], timeseries, color=ts_lines[ti][0], linestyle=ts_lines[ti][1], label=ts_labels[ti]) 
		
		plt.xlim(ts_xlim)
		plt.ylim(ts_ylim)
		plt.xticks(ts_xticks)
		plt.yticks(ts_yticks)
		plt.title(ts_title, fontsize=fontsize)
		plt.xticks(fontsize=10)
		plt.yticks(fontsize=10)
		plt.xlabel(ts_xlabel, fontsize=fontsize)
		plt.ylabel(ts_ylabel, fontsize=fontsize)
		plt.grid(color='k', linestyle='dashed', linewidth=0.5)
		if ts_insetTitle is not None:
			props = dict(boxstyle='round', facecolor='white', alpha=0.5)
			ax.text(insetTitleX, insetTitleY, ts_insetTitle, transform=ax.transAxes, fontsize=insetTitleFontsize, verticalalignment='top', bbox=props)

		plt.legend(prop={'size': ts_legendFontsize}, ncol=ts_legendNcol, loc=ts_legendLoc)
	
	#==
	
	#plt.tight_layout()
	if save:
		plt.savefig(outpath + outname)
		
	if show:
		plt.show()
			

#==
