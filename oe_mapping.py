#!/usr/bin/env python
######################################## M. Gascon. LBL 2015. ############################
# USAGE: if you have python installed (included in the PaTH)  
# open a terminal on linux, mac or windows and execute
# 
# > oe_mapping.py <filename.dat> 
#
# where filename.dat is a file from Fluorescence (Horiba)
###################################################################
 
# example of file format:  Positions (in brackets) (first line).
# wavelength and intensity (2nd line .. last line)
# Wavelength      (301,386)  (337,386)  (373,386) ...   
# ...
# 330     230     50      130     170     200     190     130     130     270 ...
# 335     160     190     190     160     140     180     160     140     150 ...
#############################################################################################  


# import libraries
import sys
import collections
import numpy as np
import re
import pandas as pd
from pandas import *

## import from glue 
from glue import qglue
from glue.core import Data, DataCollection
from glue.core.data_factories import load_data
from glue.qt.glue_application import GlueApplication
import matplotlib.pyplot as plt

#from glue.core.fitters import BaseFitter1D
#from glue.config import fit_plugin

## See number of arguments: 1st should be data and 2nd should be an image
if (len(sys.argv)>1):
	with open(sys.argv[1], 'r') as my_file:
    	       Data = my_file.readlines()
if (len(sys.argv)>2):
	with open(sys.argv[2], 'r') as my_image:
	       mapping = load_data(my_image) 

 

############################################## Read the file 
## Variables

# Data_clean will contain the data as collection of lists
Data_clean = collections.defaultdict(list)  
First_Line = 0
# Vectors with positions 
Position_x = []
Position_y = []
Wavelengths =[]

## Read the data stored in Data
for line in Data:
	line = line.rstrip()
	if not line.strip():
         continue
	else:
         Nb_Spectrum = -1
         if(First_Line == 0):             
             line = re.sub('[Wavelength\t]', '', line)
             line = re.sub('[!(]', '', line)
             line = re.sub('[!)]', '\t', line)
             Value_Column = line.split("\t") 
             for Index, Individual_Value in enumerate(Value_Column):
                if(Individual_Value.split(',')[0]):
                     Position_x.append(int(Individual_Value.split(',')[0]))
                     Position_y.append(int(Individual_Value.split(',')[1]))
             First_Line = 1
         else:
             Value_Column = line.split("\t")         
             for Index, Individual_Value in enumerate(Value_Column):
                if(Index == 0): 
                 	Wavelengths.append(int(Individual_Value))
                else: 
    			Data_clean[Index].append(float(Individual_Value))
                Nb_Spectrum = Nb_Spectrum + 1
my_file.close()
#print Data_clean

############################################### 
# indexes is a vector with numbers from 1 to the total Nb of points
indexes = [i+1 for i in xrange(len(Data_clean))]

# df is a Data Frame with columns=wavelenghts, rows=points 
df = DataFrame(index=indexes, columns=Wavelengths)

# Fill the DataFrame with rows of the Data_clean
for ix in range(len(indexes)):
	df.loc[ix+1] = Data_clean[ix+1]
#print (df)

#Position will contain the Point Number (index) and the xy positions
Positions = {'i': indexes, 'x':  Position_x , 'y':  Position_y}







#dc = DataCollection([Positions, df])

# link positional information
#dc.add_link(LinkSame(image.id['World x: RA---TAN'], catalog.id['RAJ2000']))
#dc.add_link(LinkSame(image.id['World y: DEC--TAN'], catalog.id['DEJ2000']))

#start Glue
#app = GlueApplication(dc)
#app.start()

#load_data('image.fits', factory=gridded_data)


#df.plot(df, x=Wavelengths, y=df[1])
# xyplot(df[1] ~ df[2], data = df)
#df.plot(style=['o','rx'])
# dc = DataCollection([image, catalog])

#qglue(xy=image)

#qglue(iuv=Positions, xyz=df)

#
#data_t = Data(Positions, label="positions")




#pandas_data = pd.DataFrame({'x':  Position_x, 'y':  Position_y, 'z':Data_clean})
#pandas_data = pd.DataFrame({'z':Data_clean})


#qglue(uv=Positions, xyz=pandas_data)

#fig, axes = plt.subplots(nrows=2, ncols=2)
#fig.set_figheight(6)
#fig.set_figwidth(8)
#df[0].plot(ax=axes[0,0], style='r', label='Series'); axes[0,0].set_title(0)
#df[1].plot(ax=axes[0,1]); axes[0,1].set_title(1)
#df[2].plot(ax=axes[1,0]); axes[1,0].set_title(2)
#df[3].plot(ax=axes[1,1]); axes[1,1].set_title(3)
#fig.tight_layout()

#xyplot(Data_clean[0],Data_clean[1])


'''
print "Number of Spectra: ", Nb_Spectrum

Total_gauss1_position = []
Total_gauss2_position = []
Total_gauss1_width = []
Total_gauss2_width = []
Total_gauss1_amplitude = []
Total_gauss2_amplitude = []
Total_integration = []
Gauss1_integration = []
Gauss2_integration = []

for spectrum in xrange(Nb_Spectrum+1):
  	print spectrum	
	if(spectrum > 0):
		

'''






'''
for spectrum in xrange(Nb_Spectrum+1):
    print spectrum
    if(spectrum > 0):
        gauss1  = GaussianModel(prefix='g1_')
        gauss2  = GaussianModel(prefix='g2_')
        gauss1.set_param('center',    417, min=400, max=425)
        gauss1.set_param('sigma',      15, min=10, max=20)
        gauss1.set_param('amplitude', 700000, min=10)
        gauss2.set_param('center',    425, min=420, max=430)
        gauss2.set_param('sigma',      20, min=10, max = 25)
        gauss2.set_param('amplitude', 100000, min=10)
        model = gauss1 + gauss2
        out = model.fit(np.array(Data_clean[spectrum]), x=np.array(Data_clean[0]))
        
        xx = np.linspace(np.min(Data_clean[0]), np.max(Data_clean[0]), 1000)        
        p1 = [model.params['g1_amplitude'].value, model.params['g1_center'].value, model.params['g1_sigma'].value]
        p2 = [model.params['g2_amplitude'].value, model.params['g2_center'].value, model.params['g2_sigma'].value]
        Total_gauss1_position.append(model.params['g1_center'].value)
        Total_gauss2_position.append(model.params['g2_center'].value)
        Total_gauss1_width.append(model.params['g1_fwhm'].value)
        Total_gauss2_width.append(model.params['g2_fwhm'].value)
        Total_gauss1_amplitude.append(model.params['g1_amplitude'].value)
        Total_gauss2_amplitude.append(model.params['g2_amplitude'].value)
        Total_integration.append(scipy.integrate.simps(gauss_func(xx, *p1)+gauss_func(xx, *p2)))
        Gauss1_integration.append(scipy.integrate.simps(gauss_func(xx, *p1)))
        Gauss2_integration.append(scipy.integrate.simps(gauss_func(xx, *p2)))

'''


#



#data = Data(xf=[1, 2, 3], label="first dataset")

#x_id = data.add_component( Position_x, label = 'X')
#y_id = data.add_component( Position_y, label = 'Y')
#emission_data = data.add_component(Data_clean, label = 'Emission_data')


#collection = DataCollection([data])

#print(collection)



#qglue(xyz=pandas_data)


#import matplotlib.pyplot as plt
#import scipy.interpolate
#
#
#import scipy.stats
#import scipy.optimize
#from lmfit.models import GaussianModel
#from numpy import sqrt, pi, exp


#import myfunctions

#Input file
#print sys.argv[1]


#File_to_Read = open(filename,"r")
#Data = File_to_Read.readlines()
#Data_clean = collections.defaultdict(list)
