
# Plots IGRINS spectra with models with text files with only 2 columns (BT-Settl)
# plots as many models as are in the model_path directory (i.e. file name starts with lowercase 'L')
# rotationally broadening requires evenly spaced data, so I upgrade the resultion to interpolate onto a new wavelgnth grid (newwvl) and then rotationally broaden,
# ...and THEN degrade resultion of the models to match IGRINS R
# 
# Written by Emily Lubar with many big chunks recycled from Eunkyu's original plotting code plot_igrins.py
# April 2021

import pdb
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import scipy
from scipy import stats
from astropy.io import ascii 
import glob
from os import walk
# from PyAstronomy import funcFit as fuf
from PyAstronomy import pyasl

#____________________________________________________________

b=1 #for eventual offset on plot
bb=1 #for eventaul offset between diff rotational broadenings
f = [] #lists for getting models
F = [] #lists for getting models
z = 2.26452602
xx = 0
index = 0
res_element = 3.3 #~pixels per IGRINS resolution element. This is the average res element, may nbe more or less

temp_wave_list = []
model_path = '/Users/egl637/Documents/BrownDwarfs/Emily_spectra_plotting/'
data_path = '/Users/egl637/Documents/BrownDwarfs/ABDor_L5s/'               # IGRINS data location
obs_file = 'SDCK_20201103_0048.spec_a0v.fits'           # File name of IGRINS data
# model_file = "lte012.0-3.5-0.0a+0.0.BT-Settl.spec.7.dat.txt"
low = 1.49      # Lower-end of the BT-Settl model to look, in microns. For H-band: 1.49 and for K-band: 1.85
up = 1.8        # Upper-end of the BT-Settl model to look, in microns For H-band: 1.85 and for K-band: 2.4


#____________________________________________________________

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
F.sort() #so that the models plot in descending order

#____________________________________________________________

#then for a given order, do all the rest:
for i in range(orders):
	print('step 1: iterating through orders')
	if i == 7:				#Allows you to chose one order to look at. order 9 for Aluminum region in H-band. Order 6 for CO band head in K-band
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

		plt.figure(figsize = (9, 9)) #Initialize figure before looping through models 
		
		for modelz in F:
			print('step 3: reading in model#'+str(b))
			# temp_wave, temp_flam = readBT(model, R=40000, npix=10, waverange=[low, up]) #change R for resolution: R_IGRINS = 40000
			model = ascii.read(model_path + modelz)
			temp_wave = model['col1']
			temp_flam = model["col2"]
			# print('atfirst,', type(temp_wave))

			# ### choose the indiceis of IGRINS wavelngth within the model:
			# upperlim=list(temp_wave).index(30000) ####
			# lowerlim=list(temp_wave).index(10000) ####

			# temp_flam = temp_flam[lowerlim:upperlim]
			# temp_wave = temp_wave[lowerlim:upperlim]

			temp_wave = temp_wave/1e4          # model wavelenght now in microns instead of angstroms
			
	 # Match the BT-settl model based on the IGRINS orders
			lower = np.where(np.abs(temp_wave - obs_wvl[0]) == np.min(np.abs(temp_wave - obs_wvl[0])))
			upper = np.where(np.abs(temp_wave - obs_wvl[-1]) == np.min(np.abs(temp_wave - obs_wvl[-1])))

    #stoed cut off model 
			temp_wvl_valid = temp_wave[lower[0][0]:upper[0][0]]
			temp_spec_valid = temp_flam[lower[0][0]:upper[0][0]]
	# Normalize the science and template spectra ##Move this out of loop
			obs_spec_norm = (obs_spec_valid) / np.median(obs_spec_valid)
			temp_spec_norm = (temp_spec_valid) / np.median(temp_spec_valid)

			#### make new evenly spaced wavelength grid to put everything on becuase you can't rotbraod something that is unevenly spaced:
			newwlv = np.linspace(obs_wvl[0], obs_wvl[-1], 20000) #igrins has 2000 pixles, so increase resolution but 10
			newwlvflux = np.full(20000,1)
			### interpolating model to have higher resolution to keep all information in model by going to higher resolution, and now data is evenly spaced and ready for rotbroad 
			temp_spec_norm_highR=np.interp(newwlv, temp_wvl_valid, temp_spec_norm)
			####^^this takes as arguments: the model wavelength (temp_wvl_valid), interpolates it onto the new wavelength (newwlv), and carries through the model flux (temp_spec_norm) so no flux is lost


			### now add rotational broadening to the model:
			temp_spec_norm_rotB1 = pyasl.rotBroad(newwlv, temp_spec_norm_highR, 0.75, 12/res_element)
			# temp_spec_norm_rotB2 = pyasl.rotBroad(newwlv, temp_spec_norm_highR, 0.75, 10/res_element)
			# temp_spec_norm_rotB3 = pyasl.rotBroad(newwlv, temp_spec_norm_highR, 0.75, 15/res_element)
			####^^^arguments are: wavelength to be broadened (in Angstromes), flux to be broadened, limbdarkening (chosing 0.75) and vsini in km/s (chosing 40)
			### NOTE: had newwlv/1e4 but then did some tests and the units of wavelength doesn't seem to have any impact on resulting spectra shape

			####Now that we've rotationally broadened with all the information, degrade the model to match the resolution of IGRINS data:
			temp_spec_norm_withB1=np.interp(obs_wvl, newwlv, temp_spec_norm_rotB1) ##Model wavel, observed wave, model flux gives new model flux
			# temp_spec_norm_withB2=np.interp(obs_wvl, newwlv, temp_spec_norm_rotB2) 
			# temp_spec_norm_withB3=np.interp(obs_wvl, newwlv, temp_spec_norm_rotB3)

			##plot:
			print('step 4: plotting model#'+str(b))
			plt.title('ORDER #' + str(i) +'\n'+ 'Spectra:'+data_path[36:-2] +'\n'+ ' and model(s))') #+'\n'+'RV/BC shift, res degraded, rot broad')
			# plt.plot(temp_wvl_valid+0.00068, temp_spec_norm - 1*b, label = 'Teff =' +modelz[10:11]+'00K')         # NON-ROTATIONALLY BROADENED TO COMPARE.     +5 is to offset spectra on purpose (with +0.00068 wavelength shift to account for RV)
			plt.plot(obs_wvl+0.00068,temp_spec_norm_withB1- 1*b,label = 'Teff =' +modelz[4:6]+'00K, Log(g)='+modelz[9:12]+', 12km/s broadening')  				#+5 is to offset spectra on purpose (with +0.00068 wavelength shift to account for RV)
			# plt.plot(obs_wvl+0.00068,temp_spec_norm_withB2- 2*b,label = 'Teff =' +modelz[4:6]+'00K, Log(g)='+modelz[9:12]+', 10km/s broadening')   			#+5 is to offset spectra on purpose (with +0.00068 wavelength shift to account for RV)
			# plt.plot(obs_wvl+0.00068,temp_spec_norm_withB3- 3*b,label = 'Teff =' +modelz[4:6]+'00K, Log(g)='+modelz[9:12]+', 15km/s broadening')   			#+5 is to offset spectra on purpose (with +0.00068 wavelength shift to account for RV)
			
			### to check that the interpolaton step works, plot these two (with one model at a time):
			# plt.scatter(newwlv,newwlvflux- 2*b, label = 'new evenly spaced grid',s=0.1)
			# plt.scatter(newwlv+0.00068, temp_spec_norm_highR-3*b, label = 'evenly spaced grid with fluxes',s=0.1)
			
			plt.ylabel('Normalized Flux', fontsize = 17, labelpad = 10)
			plt.xlabel('Wavelength [um]',fontsize = 17, labelpad = 10)

			b=b+1 #To add a vertical offset ot each model 
		bb=bb+1
		plt.plot(obs_wvl_valid, obs_spec_norm, label = 'IGRINS spectra')
		plt.legend(loc='upper left', fontsize = 13)
		plt.xticks(fontsize = 15)
		plt.yticks(fontsize = 15)
		# plt.ylim(-5.2,4.5) #for ~5 stacked spectra
		# plt.ylim(-7,9)

		# plt.ylim(-1,2) #for just 2 spectra
		plt.xlim(2.2925, 2.296) #for CO bandhead, order 7 in K-band
		# plt.xlim(1.673, 1.677) #for Al region, order 9 in H-band
		
		# plt.savefig('ABDOR_L5_BTsettlC_Teff_1500K_logg40_CO_zoomed.pdf')
		plt.show()




