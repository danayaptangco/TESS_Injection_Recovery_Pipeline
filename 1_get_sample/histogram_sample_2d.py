import numpy as np
from dotmap import DotMap

import matplotlib.pyplot as plt
import matplotlib.colors as colors
plt.rcParams['axes.formatter.min_exponent'] = 3

from scipy.stats import ks_2samp

class Samples():
	"""Class that contains the target, parent, and control sample"""

	def __init__(self, name, xlabel='x var [units]', 
			ylabel='y var [units]', zlabel='z var [units]'):
		''' Initialize the class and name the labels on the axes.
		'''
		self.name = name
		self.xlabel= xlabel
		self.ylabel= ylabel
		self.zlabel= zlabel

	def set_parent(self, name, xvar, yvar, extra=None, ID=None):
		''' Define the parent sample that will be drawn from
		variables x and y are used in drawing the control sample
		variable z is an optional extra variable for additional 
			verification that the control sample matches the target sample
		variable ID is an optional identifier (such as a TIC or Gaia ID)
		'''
		self.parent= DotMap()
		self.parent.name= name
		self.parent.xvar = xvar
		self.parent.yvar = yvar
		if extra is not None:
			self.parent.zvar=extra
		if ID is not None:
			self.parent.ID = ID

	def set_target(self, name, xvar, yvar, extra=None, ID=None):
		''' Define the target sample whose properties will be mimicked '''
		self.target= DotMap()
		self.target.name= name
		self.target.xvar = xvar
		self.target.yvar = yvar
		if extra is not None:
			self.target.zvar=extra
		if ID is not None:
			self.target.ID = ID

	def set_xbins(self, xmin, xmax, nx, Log=True):
		''' Define the range/grid for the x variable 
		Values outside the range will not be used in the analysis
		'''
		self.xmin= xmin
		self.xmax= xmax
		self.nx= nx
		self.xlog= Log
		if Log:
			self.xbins=np.geomspace(xmin, xmax, nx)
		else:
			self.xbins=np.linspace(xmin, xmax, nx)

	def set_ybins(self, ymin, ymax, ny, Log=True):
		''' Define the range/grid for the y variable 
		Values outside the range will not be used in the analysis
		'''
		self.ymin= ymin
		self.ymax= ymax
		self.ny= ny
		self.ylog= Log
		if Log:
			self.ybins=np.geomspace(ymin, ymax, ny)
		else:
			self.ybins=np.linspace(ymin, ymax, ny)

	def set_zbins(self, zmin, zmax, nz, Log=True):
		''' Define the range/grid for the z variable '''
		self.zmin= zmin
		self.zmax= zmax
		self.nz= nz
		self.zlog= Log
		if Log:
			self.zbins=np.geomspace(zmin, zmax, nz)
		else:
			self.zbins=np.linspace(zmin, zmax, nz)

	def plothist_x(self, Parent=True, Draw=False):
		''' plot the histogram of the two samples 
		Use the Parent and Draw boolean keywords to plot 
		either the comparison wit the parent sample,
		or the comparison with the control sample 
		'''
		self._print_in_out_x(self.target)

		if Parent:
			self._print_in_out_x(self.parent)
		if Draw:
			self._print_in_out_x(self.draw)

		if self.xlog: plt.semilogx()

		_, bins, _ = plt.hist(self.target.xvar, bins=self.xbins, 
			alpha=0.5, density=True, label=self.target.name)
		
		if Parent:
			plt.hist(self.parent.xvar, bins=self.xbins, alpha=0.5, density=True, 
				label=self.parent.name)
			
			# ks test (igoring limits)
			ks_all= ks_2samp(self.parent.xvar, self.target.xvar)
			print (ks_all)

		if Draw:
			# problem: comparison sample does not have y cut
			plt.hist(self.draw.xvar, bins=self.xbins, alpha=0.5, density=True, 
				label=self.draw.name)

		plt.xlabel(self.xlabel)

		plt.legend()
		
		#return plt

	def plothist_y(self, Parent=True, Draw=False):
		''' plot the histgoram of the two samples 
		Use the Parent and Draw boolean keywords to plot 
		either the comparison wit the parent sample,
		or the comparison with the control sample 
		'''
		self._print_in_out_y(self.target)

		if Parent:
			self._print_in_out_y(self.parent)
		if Draw:
			self._print_in_out_y(self.draw)

		if self.ylog: plt.semilogx()

		_, bins, _ = plt.hist(self.target.yvar, bins=self.ybins, 
			alpha=0.5, density=True, label=self.target.name)
		
		if Parent:
			plt.hist(self.parent.yvar, bins=self.ybins, alpha=0.5, density=True, 
				label=self.parent.name)

			# ks test (igoring limits)
			ks_all= ks_2samp(self.parent.yvar, self.target.yvar)
			print (ks_all)

		if Draw:
			plt.hist(self.draw.yvar, bins=self.ybins, alpha=0.5, density=True, 
				label=self.draw.name)

		plt.xlabel(self.ylabel)
		plt.ylabel('Count density')

		plt.legend()

	def plothist_z(self, Parent=True, Draw=False):
		''' plot the histogram of the two samples 
		Use the Parent and Draw boolean keywords to plot 
		either the comparison wit the parent sample,
		or the comparison with the control sample 
		'''

		if self.zlog: plt.semilogx()

		if 'indices' in self.target:
			target_zvar= self.target.zvar[self.target.indices]
			if Parent: parent_zvar= self.parent.zvar[self.parent.indices]
			if Draw: draw_zvar= self.draw.zvar[self.draw.indices]
		else:
			target_zvar= self.target.zvar
			if Parent: parent_zvar= self.parent.zvar
			if Draw: draw_zvar= self.draw.zvar

		_, bins, _ = plt.hist(target_zvar, bins=self.zbins, 
			alpha=0.5, density=True, label=self.target.name)
		
		if Parent:
			plt.hist(parent_zvar, bins=bins, alpha=0.5, density=True, 
				label=self.parent.name)

			# ks test
			ks_all= ks_2samp(parent_zvar, target_zvar)
			print (f'before: {ks_all}')

		if Draw:
			plt.hist(draw_zvar, bins=bins, alpha=0.5, density=True, 
				label=self.draw.name)

			# ks test
			ks_all= ks_2samp(draw_zvar, target_zvar)
			print (f'after: {ks_all}')

		plt.xlabel(self.zlabel)
		plt.ylabel('Count density')

		plt.legend()

	def hist2d_target(self):
		''' Generate and plot a 2d histogram of the target sample '''
		self.target.h2d, _, _ = np.histogram2d(self.target.xvar, self.target.yvar, 
			bins=[self.xbins,self.ybins])

		fig, ax = plt.subplots()

		if self.xlog: ax.set_xscale("log")
		if self.ylog: ax.set_yscale("log")

		plt.pcolormesh(self.xbins, self.ybins, self.target.h2d.T, 
			cmap='cubehelix_r')

		plt.title(self.target.name)
		plt.xlabel(self.xlabel)
		plt.ylabel(self.ylabel)

	def hist2d_parent(self):
		''' Generate and plot a 2d histogram of the parent sample '''
		self.parent.h2d, _, _ = np.histogram2d(self.parent.xvar, self.parent.yvar, 
			bins=[self.xbins,self.ybins])

		fig, ax = plt.subplots()

		if self.xlog: ax.set_xscale("log")
		if self.ylog: ax.set_yscale("log")

		plt.pcolormesh(self.xbins, self.ybins, self.parent.h2d.T, 
			cmap='cubehelix_r')

		plt.title(self.parent.name)
		plt.xlabel(self.xlabel)
		plt.ylabel(self.ylabel)

	def hist2d_draw(self):
		''' Generate and plot a 2d histogram of the control sample '''
		self.draw.h2d, _, _ = np.histogram2d(self.draw.xvar, self.draw.yvar, 
			bins=[self.xbins,self.ybins])

		fig, ax = plt.subplots()

		if self.xlog: ax.set_xscale("log")
		if self.ylog: ax.set_yscale("log")

		plt.pcolormesh(self.xbins, self.ybins, self.draw.h2d.T, 
			cmap='cubehelix_r')

		plt.title(self.draw.name)
		plt.xlabel(self.xlabel)
		plt.ylabel(self.ylabel)

		''' Compare matrix'''
		noff= np.sum(self.target.h2d != self.draw.h2d)
		if noff > 0:
			print (f'{noff}/{self.draw.h2d.size} are off')
		else:
			print (f'all {self.draw.h2d.size} cells exact match :)')

	def hist2d_compare(self):
		''' Compare the counts in each 2D bin and check how many are over or under/sampled'''

		nover= np.sum(self.target.h2d > self.parent.h2d)
		if nover > 0:
			print (f'{nover}/{self.parent.h2d.size} are undersampled')
		else:
			print (f'all {self.parent.h2d.size} cells oversampled :)')

		ntarget= self.target.h2d.sum(dtype=int)
		nparent= self.parent.h2d.sum(dtype=int)

		scale= 1.0 * ntarget / nparent
		ratio2d= 1.0/scale * self.target.h2d / self.parent.h2d

		fig, ax = plt.subplots()

		if self.xlog: ax.set_xscale("log")
		if self.ylog: ax.set_yscale("log")

		plt.pcolormesh(self.xbins, self.ybins, ratio2d.T, 
			cmap='seismic', norm=colors.LogNorm(vmin=0.1, vmax=10.))

		plt.colorbar()

		plt.title(f'{self.target.name} / {self.parent.name}')
		plt.xlabel(self.xlabel)
		plt.ylabel(self.ylabel)

	def make_control(self, ID= True):
		''' Draw a control sample from the parent smaple with 
		with equal counts per 2D bin as the target sample
		If the ID keyword is set, the target sample is removed 
		from the control sample such that the control sample
		and target sanmple have no overlap
		'''
		
		# generate an index for the parent sample that will be used in the draw
		self.parent.index= np.arange(self.parent.xvar.size, dtype=int)

		# generate an array that will store the draw probability of each star
		self.parent.draw_prob= np.zeros_like(self.parent.xvar)

		''' Remove target sample from parent sample '''
		if ID and 'ID' in self.parent:
			target_not_in_parent= np.isin(self.target.ID, self.parent.ID, 
				invert=True)
			print (f'missing {target_not_in_parent.sum()}/{target_not_in_parent.size}')

			parent_not_in_target= np.isin(self.parent.ID, self.target.ID, 
				invert=True)
			print (f'trimming parent sample {parent_not_in_target.sum()}/'+
					f'{parent_not_in_target.size}')

			# generate sample without target IDs
			parent_index=self.parent.index[parent_not_in_target]
			parent_xvar= self.parent.xvar[parent_not_in_target]
			parent_yvar= self.parent.yvar[parent_not_in_target]
			if 'zvar' in self.parent:
				parent_zvar= self.parent.zvar[parent_not_in_target]
			#if 'ID' in self.parent: # redundant
			parent_ID= self.parent.ID[parent_not_in_target]	

			''' check oversampling again'''
			parent_h2d, _, _ = np.histogram2d(parent_xvar, parent_yvar, 
				bins=[self.xbins,self.ybins])

			nover= np.sum(self.target.h2d > parent_h2d)
			if nover > 0:
				print (f'\n{nover}/{parent_h2d.size} cells are undersampled')
			else:
				print (f'\nall {parent_h2d.size} cells oversampled :)')


		else:
			''' Control sample may contain stars from the parent sample'''
			print (f'Using entire parent sample, {self.parent.index.size}')
			parent_index=self.parent.index
			parent_xvar= self.parent.xvar
			parent_yvar= self.parent.yvar
			if 'zvar' in self.parent:
				parent_zvar= self.parent.zvar
			if 'ID' in self.parent:
				parent_ID= self.parent.ID
		
		# generate the bin indices of each star
		ix_parent= np.digitize(parent_xvar, self.xbins, right=False)
		iy_parent= np.digitize(parent_yvar, self.ybins, right=False)

		ix_target= np.digitize(self.target.xvar, self.xbins, right=False)
		iy_target= np.digitize(self.target.yvar, self.ybins, right=False)

		i_draw= []

		# loop over all 2D bins and draw control sample from target sample
		for i in range(1,self.xbins.size):
			for j in range(1,self.ybins.size):
				# number of stars to draw in this bin
				n_target= np.sum( (i==ix_target) & (j==iy_target) ) * 10
				#assert n_target== self.target.h2d[i-1,j-1]
			
				# list of indices to draw from
				TF_parent= (i==ix_parent) & (j==iy_parent)
				n_parent= np.sum(TF_parent)
				
				# start drawing
				indices_to_draw_from= parent_index[TF_parent]
				#draw_prob= n_target/TF_parent.sum()
				# this loop can be rewritten: replace=(draw_prob>1) and if TF_parent.sum()==0
				if  TF_parent.sum() >= n_target:
					# normal sampling
					draw= np.random.choice(indices_to_draw_from, size=n_target, replace=False) #***
					
					# store draw probability of each tic
					self.parent.draw_prob[indices_to_draw_from]= n_target/TF_parent.sum()

				elif TF_parent.sum()== 0:
					# empty array, how to conserve sample size?
					print (f'  Nothing to sample at i={i},j={j}, want n={n_target}')
					draw=np.array([], dtype=int) 
				else:
					print (f'Draw with replacement {i},{j}')
					draw= np.random.choice(indices_to_draw_from, size=n_target, replace=True)
					
					# store draw probability of each tic
					self.parent.draw_prob[indices_to_draw_from]= n_target/TF_parent.sum()
				
				# store the indices from the parent sample that define the control sample
				i_draw.append(draw)

		#Store the control sample that was just drawn
		idraw=np.concatenate(i_draw)

		self.draw= DotMap()
		self.draw.name= 'MC sample'

		self.draw.xvar= self.parent.xvar[idraw]
		self.draw.yvar= self.parent.yvar[idraw]

		if 'zvar' in self.parent:
			#print ('extra variable')
			self.draw.zvar= self.parent.zvar[idraw]

		if 'ID' in self.parent:
			#print ('copying ID')
			self.draw.ID= self.parent.ID[idraw]

		#self.parent.draw_prob # was already stored in the main array

		print (f'\nsummed weights parent sample: n={self.parent.draw_prob.sum()}')			
		
		print (f'\ncontrol sample n={self.draw.xvar.size}')

	def verify_control(self):
		''' Verify that the statistical properties of the control sample 
		indeed match the target sample'''

		# indexes for sample within (x,y) grid boundaries
		self.parent.indices= \
			(self.xbins[0]<self.parent.xvar) & \
			(self.parent.xvar<=self.xbins[-1]) & \
			(self.ybins[0]<self.parent.yvar) &  \
			(self.parent.yvar<=self.ybins[-1])
		self.target.indices= \
			(self.xbins[0]<self.target.xvar) &  \
			(self.target.xvar<=self.xbins[-1]) & \
			(self.ybins[0]<self.target.yvar) &  \
			(self.target.yvar<=self.ybins[-1])
		self.draw.indices= \
			(self.xbins[0]<self.draw.xvar) &  \
			(self.draw.xvar<=self.xbins[-1]) & \
			(self.ybins[0]<self.draw.yvar) &  \
			(self.draw.yvar<=self.ybins[-1])

		# Compare target and parent sample with KS test (should not match)
		ks_before_xvar= ks_2samp(self.parent.xvar[self.parent.indices], 
				self.target.xvar[self.target.indices])
		ks_before_yvar= ks_2samp(self.parent.yvar[self.parent.indices], 
				self.target.yvar[self.target.indices])
		if 'zvar' in self.target:
			ks_before_zvar= ks_2samp(self.parent.zvar[self.parent.indices], 
					self.target.zvar[self.target.indices])

		# Compare target and control sample with KS test (should match)
		ks_after_xvar= ks_2samp(self.draw.xvar[self.draw.indices], 
				self.target.xvar[self.target.indices])
		ks_after_yvar= ks_2samp(self.draw.yvar[self.draw.indices], 
				self.target.yvar[self.target.indices])
		if 'zvar' in self.target:
			ks_after_zvar= ks_2samp(self.draw.zvar[self.draw.indices], 
					self.target.zvar[self.target.indices])

		# display the statistics 
		print ('xvar: ')
		print (f' distance: {ks_before_xvar.statistic:.3f} -> {ks_after_xvar.statistic:.3f}')
		print (f' probability: {ks_before_xvar.pvalue:.3g} -> {ks_after_xvar.pvalue:.3g}')

		print ('\nyvar: ')
		print (f' distance: {ks_before_yvar.statistic:.3f} -> {ks_after_yvar.statistic:.3f}')
		print (f' probability: {ks_before_yvar.pvalue:.3g} -> {ks_after_yvar.pvalue:.3g}')

		if 'zvar' in self.target:
			print ('\nzvar: ')
			print (f' distance: {ks_before_zvar.statistic:.3f} -> {ks_after_zvar.statistic:.3f}')
			print (f' probability: {ks_before_zvar.pvalue:.3g} -> {ks_after_zvar.pvalue:.3g}')


	def _print_in_out_x(self, sample):
		''' A helper function showing the numbet of datapoints outside of the grid '''
		nleft= np.sum(sample.xvar < self.xmin)
		nin= np.sum( (self.xmin <= sample.xvar) & (sample.xvar < self.xmax) )
		nright= np.sum(self.xmax <= sample.xvar)
		print (f'{nin}/{sample.xvar.size},  left: {nleft}, right: {nright}')

	def _print_in_out_y(self, sample):
		''' A helper function showing the numbet of datapoints outside of the grid '''
		nleft= np.sum(sample.yvar < self.ymin)
		nin= np.sum( (self.ymin <= sample.yvar) & (sample.yvar < self.ymax) )
		nright= np.sum(self.ymax <= sample.yvar)
		print (f'{nin}/{sample.yvar.size},  under: {nleft}, over: {nright}')
