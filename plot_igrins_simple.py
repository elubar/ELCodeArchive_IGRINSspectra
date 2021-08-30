
# Plots IGRINS spectra with models with text files with only 2 columns (BT-Settl)
# plots as many models as are in the model_path directory (i.e. file name starts with lowercase 'L')
# Written by Emily Lubar with many big chunks recycled from Eunkyu's original plotting code
# March 2021

import pdb
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import scipy
from scipy import stats
from astropy.io import ascii 
import glob
from os import walk

b=1 #for eventual offset on plot
f = [] #lists for getting models
F = [] #lists for getting models
temp_wave_list = []
model_path = '/Users/egl637/Documents/BrownDwarfs/Emily_spectra_plotting/'
data_path = '/Users/egl637/Documents/BrownDwarfs/ABDor_L5s/'               # IGRINS data location
obs_file = 'SDCK_20201103_0048.spec_a0v.fits'           # File name of IGRINS data
# model_file = "lte012.0-3.5-0.0a+0.0.BT-Settl.spec.7.dat.txt"
low = 1.49      # Lower-end of the BT-Settl model to look, in microns. For H-band: 1.49 and for K-band: 1.85
up = 1.8        # Upper-end of the BT-Settl model to look, in microns For H-band: 1.85 and for K-band: 2.4

#Get the data
obs_data = fits.open(data_path+obs_file)     # Read the fits file
obs_spec_full = obs_data[0].data        # Flux array
obs_wvl_full = obs_data[1].data         # Wavelength array
orders = obs_data[0].data.shape[0]      # Number of orders

# read in BT-Settl models: First get all models in directory that we want:
for (dirpath, dirnames, filenames) in walk(model_path):
	f.extend(filenames)
	break
for x in f:
	if x[0] == 'l':
		F.append(x)
#F is now a list of strings of the models in the path above and will be plotted
# print(F)
F.sort() #so that the models plot in descending order
print(F)

#then for a given order, do all the rest:
for i in range(orders):
	print('step 1: iterating through orders')
	if i == 6:				#Allows you to chose one order to look at. order 9 for Aluminum region in H-band. Order 6 for CO band head in K-band
		print('step 2: doing data and plotting data')
		index = i #chose index 9 for K band for interesting features in order 9. yes, index # = order #
		obs_spec_interest = obs_spec_full[index] #specifying order to look at (index): for flux
		obs_wvl_interest = obs_wvl_full[index] #specifying order to look at (index): for wavelength

    # Only take non-NaN's
		valid = np.where(np.isfinite(obs_spec_interest)) #no more NaNs
		obs_spec_valid = obs_spec_interest[valid]
		obs_wvl_valid = obs_wvl_interest[valid]

    # Clip the outliers (7 sigma below and 5 sigma above the average)
		std = np.std(obs_spec_valid)
		mean = np.mean(obs_spec_valid)
		obs_spec_clipped, low, high = scipy.stats.sigmaclip(obs_spec_valid, low = 7.0, high = 5.0)
		array = np.where(np.logical_and(obs_spec_valid > low, obs_spec_valid < high))
		obs_wvl_clipped = obs_wvl_valid[array]

		obs_spec = obs_spec_clipped
		obs_wvl = obs_wvl_clipped

		plt.figure(figsize = (15, 7)) #Initialize figure before looping through models 
		
		for modelz in F:
			print('step 3: reading in model#'+str(b))
			# temp_wave, temp_flam = readBT(model, R=40000, npix=10, waverange=[low, up]) #change R for resolution: R_IGRINS = 40000
			model = ascii.read(model_path + modelz)
			temp_wave = model['col1']
			temp_flam = model["col2"]

# print(wvl)
			temp_wave = temp_wave/1e4          # Angstroms to Microns

	 # Match the BT-settl model based on the IGRINS orders
			lower = np.where(np.abs(temp_wave - obs_wvl[0]) == np.min(np.abs(temp_wave - obs_wvl[0])))
			upper = np.where(np.abs(temp_wave - obs_wvl[-1]) == np.min(np.abs(temp_wave - obs_wvl[-1])))

    #stoed cut off model 
			temp_wvl_valid = temp_wave[lower[0][0]:upper[0][0]]
			temp_spec_valid = temp_flam[lower[0][0]:upper[0][0]]
	# Normalize the science and template spectra ##Move this out of loop
			obs_spec_norm = (obs_spec_valid) / np.median(obs_spec_valid)
			temp_spec_norm = (temp_spec_valid) / np.median(temp_spec_valid)

			print('step 4: plotting model#'+str(b))
			plt.title('ORDER #' + str(i) +'\n'+ 'Spectra:'+data_path[36:-2] +'\n'+ ' and model(s)')
			plt.plot(temp_wvl_valid, temp_spec_norm - 1*b, label = 'Model:' +modelz[0:28])         # +5 is to offset spectra on purpose
			plt.xlabel('Wavelength ($\mu m$)', fontsize = 17, labelpad = 10)
			plt.ylabel('Normalized Flux', fontsize = 17, labelpad = 10)

			b=b+1 #To add a vertical offset ot each model 

		plt.plot(obs_wvl_valid, obs_spec_norm, label = 'IGRINS spectra')
		plt.legend(loc=0, fontsize = 13)
		plt.xticks(fontsize = 15)
		plt.yticks(fontsize = 15)
		plt.ylim(-5.2,4.5)
		plt.xlim(2.294, 2.3)
		# plt.xlim(1.673, 1.677) #for Al region, order 9 in H-band
		plt.savefig('ABDOR_L5_BTsettlC_logg55_Teff_1200K_1600K_CO.pdf')
		plt.show()
		


