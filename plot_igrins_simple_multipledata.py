
# FINISHED this will plot multiple files (not necessariuly with models, but just to have multiple objects plotted against 4each other to show temp grad) 
# but you have to specify which ones by hand which should be sufficient and allow for manueving 

# written to plot data so I can match spectral types with file nanes inthe gemini observations

# other notes:
# Plots IGRINS spectra with models with text files with only 2 columns (BT-Settl)
# plots as many models as are in the model_path directory (i.e. file name starts with lowercase 'L')
# Written by Emily Lubar with many recycled pieces from Eunkyu's original plotting code
# April 26 2021

import pdb
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import scipy
from scipy import stats
from astropy.io import ascii 
import glob
from os import walk

b=0 #for eventual offset on plot and then appropirated for other things
bb = 0#for eventual offset on plot
f = [] #lists for getting models
F = [] #lists for getting models
temp_wave_list = []

fluxes = []
wavels = []
#__________________________________________________

# recipefile= 20210101
# specttype = ['G2','K2','K4','K4']
# obs_file1 = 'SDCH_20210101_0235.spec_a0v.fits'           #0101 observing: ordered
# obs_file2 = 'SDCH_20210101_0212.spec_a0v.fits'
# obs_file3 = 'SDCH_20210101_0174.spec_a0v.fits'
# obs_file4 = 'SDCH_20210101_0170.spec_a0v.fits'

# recipefile= 20210102
# specttype = ['G8','M9','M5','M4.3','M4','M3.5','M3.5','M3','L5']
# obs_file1 = 'SDCH_20210102_0041.spec_a0v.fits'           # 0102 observing: ________
# obs_file2 = 'SDCH_20210102_0071.spec_a0v.fits'
# obs_file3 = 'SDCH_20210102_0063.spec_a0v.fits'
# obs_file4 = 'SDCH_20210102_0097.spec_a0v.fits'
# obs_file5 = 'SDCH_20210102_0130.spec_a0v.fits'
# obs_file6 = 'SDCH_20210102_0055.spec_a0v.fits'
# obs_file7 = 'SDCH_20210102_0081.spec_a0v.fits'
# obs_file8 = 'SDCH_20210102_0047.spec_a0v.fits'
# obs_file9 = 'SDCH_20210102_0110.spec_a0v.fits' ##this L5 is a little weird and noisy


# recipefile= 20210103
# obs_file1 = 'SDCH_20210102_0041.spec_a0v.fits'           # 0102 observing: ordered
# obs_file2 = 'SDCH_20210102_0047.spec_a0v.fits'

# recipefile= 20210104
# obs_file1 = 'SDCH_20210104_0021.spec_a0v.fits'           # 0102 observing: ordered

# recipefile= 20210110
# obs_file1 = 'SDCH_20210110_0021.spec_a0v.fits'           # 0102 observing: ________
# obs_file2 = 'SDCH_20210110_0025.spec_a0v.fits'   

# recipefile= 20210111
# obs_file1 = 'SDCH_20210104_0021.spec_a0v.fits'           # 0102 observing: ordered

# recipefile= 20210115
# obs_file1 = 'SDCH_20210115_0021.spec_a0v.fits'  


model_path = '/Users/egl637/Documents/BrownDwarfs/Emily_spectra_plotting/'
# data_path = '/Users/egl637/Documents/BrownDwarfs/ABDor_L5s/'               # IGRINS data location
# data_path = '/Users/egl637/Documents/BrownDwarfs/GS-2020B-Q-319_Lubar/'+str(recipefile)+'/reduced/'               # IGRINS data location
# obsfileS = [obs_file1,obs_file2,obs_file3,obs_file4,obs_file5,obs_file6,obs_file7,obs_file8,obs_file9]
#__________________________________________________


low = 1.49      # Lower-end of the BT-Settl model to look, in microns. For H-band: 1.49 and for K-band: 1.85
up = 1.8        # Upper-end of the BT-Settl model to look, in microns For H-band: 1.85 and for K-band: 2.4

choseband = 'H'
chosenorder = 9
#OH doublet is order 17 in H band
#Ca II is order 10 in H band
#Al is order 9 in H

obsfileS = []
data_path = '/Users/egl637/Documents/BrownDwarfs/GS-2020B-Q-319_Lubar/All_ABDorsample'+choseband+'/'

for (dirpath, dirnames, filenames) in walk(data_path):
	f.extend(filenames)
	break
for x in f:
	if x[0] == 'S':
		obsfileS.append(x)
#F is now a list of strings of the models in the path above and will be plotted
# F.sort() #so that the models plot in descending order
print(obsfileS)
print('length1',len(obsfileS))

###The following two lists are made with the "AB Dor GS-2020B-Q-319 observation mapping: files to names" google sheet:
# obsfileS = [obsfileS[12],obsfileS[6],obsfileS[3],obsfileS[10],obsfileS[9],obsfileS[1],obsfileS[15],obsfileS[0],obsfileS[7],obsfileS[4],obsfileS[2],obsfileS[13],obsfileS[17],obsfileS[8],obsfileS[14],obsfileS[11],obsfileS[5],obsfileS[16]]

obsfileS = [obsfileS[16],obsfileS[5],obsfileS[11],obsfileS[14],obsfileS[8],obsfileS[17],obsfileS[13],obsfileS[2],obsfileS[4],obsfileS[7],obsfileS[0],obsfileS[15],obsfileS[1],obsfileS[9],obsfileS[10],obsfileS[3],obsfileS[6],obsfileS[12]]
# specttype = ['G2', 'G8', 'K2','K3','K4','K4','M0.5','M3','M3','M3','M3.5','M3.5','M4','M4.3', 'M4.5', 'M5', 'M9','L5']
specttype = ['L5','M9','M5','M4.5','M4.3','M4','M3.5','M3.5','M3','M3','M3','M0.5','K4','K4','K3','K2','G8', 'G2']

#__________________________________________________


# read in BT-Settl models: First get all models in directory that we want:
for (dirpath, dirnames, filenames) in walk(model_path):
	f.extend(filenames)
	break
for x in f:
	if x[0] == 'l':
		F.append(x)
#F is now a list of strings of the models in the path above and will be plotted
F.sort() #so that the models plot in descending order
print('length2',len(F))


#Get the data
for theobjectfile in obsfileS:
	obs_data = fits.open(data_path+theobjectfile)     # Read the fits file
	obs_spec_full = obs_data[0].data        # Flux array
	obs_wvl_full = obs_data[1].data         # Wavelength array
	orders = obs_data[0].data.shape[0]      # Number of orders



#then for a given order in the data:
	for i in range(orders):
		print('step 1: iterating through orders')
		if i == chosenorder:				#Allows you to chose one order to look at. order 9 for Aluminum region in H-band. Order 6 for CO band head in K-band
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

			# for modelz in F:
			print('step 3: reading in model#'+str(b))

			model = ascii.read(model_path + F[0])
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

			fluxes.append(obs_spec_norm)
			wavels.append(obs_wvl_valid)
			print('just finished orders loop, '+str(b))

# Now plot all the wavelengths and fluxes of each object

plt.figure(figsize = (7, 7)) #Initialize figure before looping through data 


print('number of objects being plotted:', len(fluxes))
print(b)

# plt.axvline(1.649,-2,6, label = 'Ca II line @ 1.649 um')
# plt.axvline(1.651,-2,6, label = 'Ca II line @ 1.651 um')
# plt.axvline(1.56265,-2,6, label = 'OH doublet')


for r,e in zip(wavels,fluxes):

	plt.title('ORDER #' + str(chosenorder) +'\n'+ 'Spectra:'+data_path[36:] +'\n'+ ' and model(s)')
		# plt.plot(temp_wvl_valid, temp_spec_norm - 1*b, label = 'Model:' +F[0][0:28])         # +5 is to offset spectra on purpose
	plt.xlabel('Wavelength ($\mu m$)', fontsize = 17, labelpad = 10)
	plt.ylabel('Normalized Flux', fontsize = 17, labelpad = 10)
	print('before', b)
	# plt.plot(r, e+b, label = 'IGRINS spectra for object:{}'.format(obsfileS[b][14:18])+' spectral type:'+specttype[b])
	plt.plot(r, e+bb, label = specttype[b])

	## to reverse the labels in the ledgend so they match the order of the plot:
	current_handles, current_labels = plt.gca().get_legend_handles_labels()
		# sort or reorder the labels and handles
	reversed_handles = list(reversed(current_handles))
	reversed_labels = list(reversed(current_labels))
	# call plt.legend() with the new values
	plt.legend(reversed_handles,reversed_labels,loc = 'center left',bbox_to_anchor=(1, 0.5), fontsize = 13, title='spectral type')
	#Then to keep the legend on the screen for a given fig size:
	plt.tight_layout(rect=[0,0,0.75,1.0])

	#futher adjustments
	plt.xticks(fontsize = 15)
	plt.yticks(fontsize = 15)
	plt.ylim(0,20)
	# plt.xlim(2.294, 2.3)
	# plt.xlim(1.671, 1.679) #for Al region, order 9 in H-band
	# plt.xlim(1.560, 1.566) #for Al region, order 9 in H-band

	# plt.savefig('ABDOR_L5_BTsettlC_logg55_Teff_1200K_1600K_CO.pdf')

	b=b+1 #To add a vertical offset ot each 
	bb = bb+0.5
# plt.plot(obs_wvl+0.00068,temp_spec_norm- 0.33 ,label = 'Teff =' +F[0][4:6]+'00K, Log(g)='+F[0][9:12])

plt.show()

print('after', b)
		

