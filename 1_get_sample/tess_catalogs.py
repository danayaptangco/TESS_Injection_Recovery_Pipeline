import numpy as np
import pandas as pd
import os, sys #, glob
from dotmap import DotMap

from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
# print (NasaExoplanetArchive.TAP_TABLES)

def astroquery_TOIs(date = None, Save_to_file=False):
	""" Fetches the TOIs from the exoplanet archive.
	Note: only saves a couple stellar parameters
	Note: discards all NaN entries as well
	if you set the date, it will save the catalogue to a file for future use
	if you set the Save_to_file keyword, it will re-generate the catalogue
	"""

	toi= DotMap() 

	if date is not None:
		fname= f'TOI_{date}.npz'
		if os.path.isfile(fname):
			print (f'Read from file: {fname}')
			npz= np.load(fname)
			toi.Rst= npz['Rst']
			toi.Teff= npz['Teff']
			toi.d= npz['d']
			toi.Tmag= npz['Tmag']
			toi.tic= npz['tic']
			npz.close()
		else:
			Save_to_file= True



	if (date is None) or Save_to_file:
		TOI= NasaExoplanetArchive.query_criteria(table='toi')
		# columns: TOI.colnames

		# get rid of astropy masked array bullshit
		toi.Rst= TOI['st_rad'].filled(np.nan).value
		toi.Teff= TOI['st_teff'].filled(np.nan).value
		#logg= TOI['st_logg']
		#L= TOI['st_']
		toi.d= TOI['st_dist'].filled(np.nan).value
		toi.Tmag= TOI['st_tmag'].filled(np.nan).value

		toi.tic= TOI['tid'].filled().value

		nope= np.isfinite(toi.Rst) & np.isfinite(toi.Teff) & np.isfinite(toi.d) & np.isfinite(toi.Tmag)

		print (nope.size, nope.sum())
		toi.Rst= toi.Rst[nope]
		toi.Teff= toi.Teff[nope]
		toi.d= toi.d[nope]
		toi.Tmag= toi.Tmag[nope]
		toi.tic= toi.tic[nope]

		# unique?
		u, indices = np.unique(toi.tic, return_index=True)
		if toi.tic.size  > u.size:
			print (f'duplicate tics {u.size}/{toi.tic.size } unique')
			toi.Rst= toi.Rst[indices]
			toi.Teff= toi.Teff[indices]
			toi.d= toi.d[indices]
			toi.Tmag= toi.Tmag[indices]
			toi.tic= toi.tic[indices]			

		if Save_to_file:
			print (f'Save to file: {fname}')
			np.savez_compressed(fname, Rst=toi.Rst, Teff=toi.Teff, 
				d=toi.d, Tmag=toi.Tmag, tic=toi.tic)

	return toi

def CTL(path= './', Gaia=False, DropNaN=True, date = None, Save_to_file=False):
	''' Read the CTL using astropy or from file
	Note: only saves a couple stellar parameters
	Note: discards all NaN entries dith DropNaN
	Note: The Gaia ID is optional because it causes some issues with strings that are NaNs
	if you set the date, it will save the catalogue to a file for future use
	if you set the Save_to_file keyword, it will re-generate the catalogue
	'''

	''' previously used/saved version'''
	if date is not None:
		fname= f'CTL_{date}.npz'
		if os.path.isfile(fname):
			print (f'Read from file: {fname}')
			npz= np.load(fname)

			ctl=DotMap() 
			ctl.Rst= npz['Rst']
			ctl.Teff= npz['Teff']
			ctl.Mst=npz['Mst']
			ctl.Lst=npz['Lst']
			ctl.d= npz['d']
			ctl.Tmag= npz['Tmag']
			ctl.tic= npz['tic']
			if 'gaia' in npz:
				ctl.gaia=npz['gaia']
			npz.close()
		else:
			Save_to_file= True

	if (date is None) or Save_to_file:
		if Gaia:
			# Force the data type instead
			dtype= {'[ID]:Integer':int,
				'[GAIA]:String':str,
				'[Tmag]:Float':float,
				'[GAIAmag]:Float':float,
				'[Teff]:Float':float,
				'[rad]:Float':float,
				'[mass]:Float':float,
				'[lum]:Float':float,
				'[d]:Float':float}
			columns=list(dtype)

			df= pd.read_csv(path+"TESS/CTL/CTL_with_header.csv", usecols=columns, dtype=dtype)
		else:
			columns= ['[ID]:Integer','[GAIA]:String',
				'[Tmag]:Float','[GAIAmag]:Float',
				'[Teff]:Float','[rad]:Float','[mass]:Float',
				'[lum]:Float','[d]:Float']
			df= pd.read_csv(path+"TESS/CTL/CTL_with_header.csv", usecols=columns)

		len_df=len(df)
		print (f'CTL has {len_df} entries')

		print ('\nNaNs: ')
		print (df.isna().sum())

		if DropNaN:
			nancolumns= ['[ID]:Integer','[GAIA]:String',
				'[rad]:Float', '[Teff]:Float','[d]:Float', '[Tmag]:Float']
			df.dropna(subset= nancolumns, inplace=True)
			print (f'Dropped {len_df-len(df)} entries')

		ctl=DotMap() 

		# get rid of astropy masked array bullshit
		ctl.Rst= df['[rad]:Float'].to_numpy()
		ctl.Teff= df['[Teff]:Float'].to_numpy()
		ctl.Mst= df['[mass]:Float'].to_numpy()
		ctl.Lst= df['[lum]:Float'].to_numpy()
		ctl.d= df['[d]:Float'].to_numpy()
		ctl.Tmag= df['[Tmag]:Float'].to_numpy()
		ctl.tic= df['[ID]:Integer'].to_numpy()
		if Gaia:
			if DropNaN:
				ctl.gaia= df['[GAIA]:String'].to_numpy(dtype=int)
			else:
				print ('Note: Gaia ID is string because NaNs')
				ctl.gaia= df['[GAIA]:String'].to_numpy()

		if Save_to_file:
			print (f'Save to file: {fname}')
			np.savez_compressed(fname, **ctl)
			#np.savez_compressed(fname, Rst=ctl.Rst, Teff=ctl.Teff, 
			#	Mst=ctl.Mst, Lst=ctl.Lst, d=ctl.d, Tmag=ctl.Tmag, tic=ctl.tic)

	return ctl