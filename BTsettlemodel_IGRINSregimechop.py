
# Takes (BT-Settl) models with text files with only 2 columns and makes new .txt files with a subset of the data that matches the IGRINS wavelengths 
# Otherwise it takes forever to process them when I run plot_igrins_simple_res_rot_v2.py

# it works!

# Written by Emily Lubar with some help from Sam
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
from PyAstronomy import pyasl

#____________________________________________________________

f = [] #lists for getting models
F = [] #lists for getting models
dd=0

temp_wave_list = []
model_path = '/Users/egl637/Documents/BrownDwarfs/Emily_spectra_plotting/'
# model_file = "lte012.0-3.5-0.0a+0.0.BT-Settl.spec.7.dat.txt"
low = 1.49      # Lower-end of the BT-Settl model to look, in microns. For H-band: 1.49 and for K-band: 1.85
up = 1.8        # Upper-end of the BT-Settl model to look, in microns For H-band: 1.85 and for K-band: 2.4

#____________________________________________________________

# read in BT-Settl models: First get all models in directory that we want:
for (dirpath, dirnames, filenames) in walk(model_path):
	f.extend(filenames)
	break
for x in f:
	if x[0] == 'l':
		F.append(x)

#F is now a list of strings of the models in the path above and will be plotted
F.sort() #so that the models plot in descending order

for i in F:
	print(i)
# sanity check to see what I'm processing

for modelz in F:
	dd = dd+1
	print('model #'+str(dd)) #to keep track if I'm running more than a few cause it takes many minutes
	model = ascii.read(model_path + modelz)
	temp_wave = model['col1']
	temp_flam = model["col2"]

	### choose the indiceis of IGRINS wavelngth within the model:
	upperlim=list(temp_wave).index(24000) ####
	lowerlim=list(temp_wave).index(15000) ####

	# print(upperlim)
	# print(lowerlim)

	###Santiy check:
	print('flux_i', len(temp_flam))
	print('wave_i', len(temp_wave))

	temp_flam = temp_flam[lowerlim:upperlim]
	temp_wave = temp_wave[lowerlim:upperlim]
	 
	###Santiy check #2:
	print('flux_f', len(temp_flam))
	print('wave_f', len(temp_wave))

	##The name of the new file will be the same as original but with an added _IGRINSwavelens.txt:
	newfilename = modelz[:-4]+'_IGRINSwavelens.txt'

	##Write the two columns of IGRINS range wavelengths to a new text file with the name specified above:
	ascii.write([temp_wave,temp_flam],newfilename) #,overwrite=True)

#Done! :)






