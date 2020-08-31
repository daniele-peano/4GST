""" 
Script to compute LAI phenology, starting from the 
method by Toratarolo et al. 2013

contact:
Daniele Peano: daniele.peano@cmcc.it
"""

# import required libraries

from netCDF4 import Dataset
import numpy as np
from datetime import datetime
from mpl_toolkits.basemap import Basemap, cm
import matplotlib.pyplot as plt
from scipy import stats
import sys
import getopt

# read option given by user

argv = sys.argv[1:]
if len(argv) == 0:
	print '*** MISSING OPTION ***'
	print 'Usage: 4GST_phenology.py -m <model> -y <year>'
try:
	opts, args = getopt.getopt(argv,"hm:y:",["model=","year="])

except getopt.GetoptError:
	print 'Usage: 4GST_phenology.py -m <model> -y <year>'
	sys.exit(2)

for opt, arg in opts:
	if opt == '-h':
		print 'Usage: 4GST_phenology.py -m <model> -y <year>'
		print '-m, --model: model to be used'
		print '-y, --year: year to be analized'
		sys.exit()
	elif opt in ("-m", "--model"):
		mdl = arg
	elif opt in ("-y", "--year"):
		anno = arg

# define the file to be used

if mdl == "LAI3g" and anno == "1982-2011":
	filename = './LAI3g/LAI3g_1982-2011_global1deg.nc' 
elif mdl == "LAI3g":
	filename = '../../LAI3G/final/LAI3g_'+anno+'_global1deg.nc' 
elif mdl == "CLM45" and anno == "1982-2011":
	filename = './CLM45/CLM45.CMCC_S3_1deg_1982-2011.nc' 
elif mdl == "CLM45" and anno == "1996-2005":
	filename = './CLM45/CLM45.CMCC_S3_1deg_1996-2005.nc' 
elif mdl == "CLM45":
	filename = '../../MODELS/CLM45_1deg/CLM45.CMCC_S3_1deg_'+anno+'.nc' 
elif mdl == "LAI3gNH" and anno == "1982-2011":
	filename = './LAI3g/LAI3g_1982-2011_NH1deg.nc' 
elif mdl == "CLM45s3b":
	filename = './CLM45s3b/CLM45.CMCC_S3b_15d_1996-2005.nc'

if mdl == "LAI3g" and anno == "1982-2011":
	outname1 = './LAI3g/Phenol_LAI3g_1982-2011_Global_1deg_4GST.nc'
elif mdl == "LAI3g":
	outname1 = './LAI3g/Phenol_LAI3g_'+anno+'_Global_1deg_4GST.nc'
elif mdl == "CLM45" and anno == "1982-2011":
	outname1 = './CLM45/Phenol_CLM45_1982-2011_Global_1deg_4GST.nc'
elif mdl == "CLM45" and anno == "1996-2005":
	outname1 = './CLM45/Phenol_CLM45_1996-2005_Global_1deg_4GST.nc'
elif mdl == "CLM45":
	outname1 = './CLM45/Phenol_CLM45_'+anno+'_Global_1deg_4GST.nc'
elif mdl == "LAI3gNH" and anno == "1982-2011":
	outname1 = './LAI3g/Phenol_LAI3g_1982-2011_NH_1deg_4GST.nc'
elif mdl == "CLM45s3b":
	outname1 = './CLM45s3b/Phenol_CLM45_S3b_1996-2005_4GST.nc'
	

#*****************************************
# READ NETCDF FILE
#*****************************************

print 'Reading:  ', filename

dataset = Dataset(filename, mode='r')

time = dataset.variables['time'][:] 
if mdl == "LAI3g" or mdl == "LAI3gNH":
	LAIr = dataset.variables['LAI'][:] # LAI3g case
	lons = dataset.variables['longitude'][:]    
	lats = dataset.variables['latitude'][:]    
elif mdl == "CLM45":
	LAIr = dataset.variables['dlai'][:] # CLM45 case
	lons = dataset.variables['longitude'][:]
	lats = dataset.variables['latitude'][:]
elif mdl == "CLM45s3b":
	LAIr = dataset.variables['TLAI'][:]
	lons = dataset.variables['lon'][:]
        lats = dataset.variables['lat'][:]

print len(lons)
print len(lats)

LAI1 = np.zeros(shape=(len(time),len(lats),len(lons)))
for nlon in range(len(lons)):
        for nlat in range(len(lats)):
                for ntim in range(len(time)):
                        if LAIr[ntim,nlat,nlon] >= 20:
                                LAI1[ntim,nlat,nlon] = np.nan
                        else:
                                LAI1[ntim,nlat,nlon] = LAIr[ntim,nlat,nlon]

dataset.close()

#*****************************************
# COMPUTE ONSET/OFFSET
#*****************************************

# ******* compute time minimum ******
LAImin = np.nanmin(LAI1, axis=0)

# ******* compute time maximum ******
LAImax = np.nanmax(LAI1, axis=0)

# ****** compute annual mean ******
LAImean = np.nanmean(LAI1, axis=0)

# ****** compute Seasonal Amplitude *****
SeasAmpl = LAImax - LAImin

# ******** Define the Onset, Offset, and seasonal type matrices ****

Onsets = np.full((2,len(lats),len(lons)), np.nan)
Offsets = np.full((2,len(lats),len(lons)), np.nan)
SeasType = np.full((len(lats),len(lons)), np.nan)
SeasClss = np.full((len(lats),len(lons)), np.nan)
SeasClss1 = np.full((len(lats),len(lons)), np.nan)

# ******** Compute the "Type" of LAI seasonal cycle ****
# Type 1) evergreen phenology case          
#           __       __
#      ____/  \   __/  \_____
#              \_/
#
# Type 2) basic case with one growing season
#             _____
#            /     \
#           /       \
#      ____/         \_____
#
# Type 3.1) two growing season, in the same year
#             __
#            /  \    __
#           /    \__/  \
#      ____/            \__
#
# Type 3.4) two growing season, in two years
#             __
#      _     /  \     __
#       \   /    \   /
#        \_/      \_/    
#
#
# This method divides the annual LAI timeseries in 4 pieces
# each of 3 months, with one overlapping step.
# Then, the linear regression is computed for each of the 4 pieces
# and a Seasonal Type is assigned depending on the directions 
# (positive or negative linear regression gradient) of the 
# four linear regressions 

# ********** Divide the time axis in 4 pieces *******

time1 = np.zeros(shape=(7))
time2 = np.zeros(shape=(7))
time3 = np.zeros(shape=(7))
time4 = np.zeros(shape=(7))

for i in range(7):
	time1[i] = time[i]
	time2[i] = time1[i]+6
	time3[i] = time1[i]+12
	time4[i] = time1[i]+18
	
# ********** Divide the LAI timeseries in 4 pieces *******

for nlon in range(len(lons)):
	for nlat in range(len(lats)):
		if LAImean[nlat,nlon] > 0:
			day = [15,30,45,60,75,90,105,120,135,150,165,180,195,210,225,240,255,270,285,300,315,330,345,360]
			LAI1dr = LAI1[:,nlat,nlon]
			knan = 0
		# identify LAI maximum 
			for t in range(len(LAI1dr)):
				if np.isnan(LAI1dr[t]):
					LAI1dr[t] = 0
					knan = knan + 1
				if LAI1dr[t] == LAImax[nlat,nlon]:
			                midPos = t
			LAI1d = np.zeros(shape=(len(LAI1dr)))
			days= np.zeros(shape=(len(day)))
		# shift of time series to place maximum LAI in the time series center
			if midPos == 11:
			        LAI1d = LAI1dr
			        days = day
			elif midPos > 11:
			        dMP = midPos - 11
			        for l in range(len(LAI1dr)-dMP):
			                LAI1d[l]=LAI1dr[l+dMP]
			                days[l]=day[l+dMP]
			        for l in range(dMP):
			                LAI1d[len(LAI1dr)-dMP+l] = LAI1dr[l]
			                days[len(day)-dMP+l] = day[l]
			else:
			        dMP = 11-midPos
			        for l in range(dMP):
			                LAI1d[l]=LAI1dr[len(LAI1dr)-dMP+l]
			                days[l]=day[len(days)-dMP+l]
			        for l in range(len(LAI1dr)-dMP):
			                LAI1d[dMP+l] = LAI1dr[l]
			                days[dMP+l] = day[l]

			LAI1d_1 = np.zeros(shape=(7))
			LAI1d_2 = np.zeros(shape=(7))
			LAI1d_3 = np.zeros(shape=(7))
			LAI1d_4 = np.zeros(shape=(7))
			for i1 in range(7):
			        i2 = i1+6
			        i3 = i1+12
			        i4 = i1+18
				LAI1d_1[i1] = LAI1d[i1]
				LAI1d_2[i1] = LAI1d[i2]
				LAI1d_3[i1] = LAI1d[i3]
				if i1 == 6:
					LAI1d_4[i1] = LAI1d[0]
				else:
					LAI1d_4[i1] = LAI1d[i4]
#			********** set to zero the nan points  *******
			for i in range(7):
			        if np.isnan(LAI1d_1[i]):
			                LAI1d_1[i] = 0
			        if np.isnan(LAI1d_2[i]):
			                LAI1d_2[i] = 0       
			        if np.isnan(LAI1d_3[i]):
			                LAI1d_3[i] = 0
			        if np.isnan(LAI1d_4[i]):
			                LAI1d_4[i] = 0

#			********** compute the 4 linear regressions **
			linreg1 = stats.linregress(time1,LAI1d_1)
			linreg2 = stats.linregress(time2,LAI1d_2)
			linreg3 = stats.linregress(time3,LAI1d_3)
			linreg4 = stats.linregress(time4,LAI1d_4)

# ********** Setting the tipe depending on the linear regression gradient sign
# Type 3.1) / \ / \; 
# Type 3.4) \ / \ /;
# Type 2) - ... ; ... > half points; else
# Type 1) evergreen areas:
# regions where seasonal LAI amplitude is lower than 25 percentage of LAI minimum

			if SeasAmpl[nlat,nlon] < 0.25*LAImean[nlat,nlon]:
				SeasType[nlat,nlon] = 1
			else:
				if linreg1[0] == 0.0:
					SeasType[nlat,nlon] = 2
				elif knan >= len(LAI1dr)/2:
					SeasType[nlat,nlon] = 2
				elif linreg1[0] >= 0 and linreg2[0] <= 0 and linreg3[0] >= 0 and linreg4[0] <= 0:
					SeasType[nlat,nlon] = 3
					SeasClss[nlat,nlon] = 1
				elif linreg1[0] <= 0 and linreg2[0] >= 0 and linreg3[0] <= 0 and linreg4[0] >= 0:
					SeasType[nlat,nlon] = 3
					SeasClss[nlat,nlon] = 4
				else:
					SeasType[nlat,nlon] = 2
				
# ******* compute the Onset and Offset in each gridpoint *****

			Ab = LAI1d
			for t in range (len(Ab)):
				if np.isnan(Ab[t]):
					Ab[t] = 0
			if SeasType[nlat,nlon] == 2:
				mx = np.max(Ab)
				mn = np.min(Ab)
				thrs = mn + 0.2 * ( mx - mn )
				if mx < 0.01:
					Onset = 0
					Offset = 0
				else:
					tr = np.where(Ab >= thrs)
					trar = tr[0]
					if len(trar) == 24:
						Onset = 0
						Offset = 0
					else:
		# operations to smooth the line and avoid one-step perturbation to disturb onset and offset detection
						kont = np.zeros(shape=(len(trar)))
						for i in range(len(trar)-1):
	        					if trar[i+1]-trar[i] <= 2: 
	               						kont[i] = 1
	        					elif trar[i+1]-trar[i] <= 2:
	                					kont[i] = 1
	       						else:
	                					kont[i] = 0
						if trar[len(trar)-1]-trar[len(trar)-2] <= 2:
						        kont[len(trar)-1] = 1
						else:
	        					kont[len(trar)-1] = 0
	
						kont2 = 0
						for i in range(len(kont)):
	        					if kont[i] == 0:
	                					kont2 = kont2 + 1
	                					break
	        					else:
	                					kont2 = kont2 + 1
	
						if kont2 == len(trar):
	        					trar2 = np.zeros(shape=(kont2))
	        					for i in range(kont2):
	                					trar2[i] = trar[i]
						elif kont[0] == 0:
	        					trar2 = np.zeros(shape=(len(kont)-1))
	        					for i in range(len(kont)-1):
	                					trar2[i] = trar[i+1]
						else:
	        					valbc2=0
	        					valac2=0
	        					for i in range(len(kont)):
	                					if i <= kont2-1:
	                        					valbc2 = valbc2+1
	                					else:
	                        					valac2 = valac2+1
	        					if valbc2 > valac2:
	                					trar2 = np.zeros(shape=(kont2))
	                					for i in range(len(trar2)):
	                        					trar2[i] = trar[i]
	        					else:
	                					trar2 = np.zeros(shape=(valac2))
	                					for i in range(len(trar2)):
	                        					trar2[i] = trar[i+kont2]
	
						Onset = trar2[0]
						Offset = trar2[len(trar2)-1]

						Onsets[0,nlat,nlon] = days[int(Onset)]
                                                Offsets[0,nlat,nlon] = days[int(Offset)]

						if lats[nlat] >= 0:
                                                        if days[11]>=120 and days[11]<=270:
                                                                SeasClss1[nlat,nlon] = 1
                                                        else:
                                                                SeasClss1[nlat,nlon] = 2
                                                else:
                                                        if days[11]>=120 and days[11]<=270:
                                                                SeasClss1[nlat,nlon] = 2
                                                        else:
                                                                SeasClss1[nlat,nlon] = 1
						
	
			elif SeasType[nlat,nlon] == 3 and SeasClss[nlat,nlon] == 1:

                        # determine the central minimum that divides the two cycles
                                LAI_cnt = Ab[7:18]
                                C31min = np.min(LAI_cnt)
                                for i in range(len(Ab)):
                                        if Ab[i] == C31min and i >= 7 and i < 18:
                                                brk_pos = i

			# analyse the first cycle as a shorter type 1 cycle
                                LAIc1 = Ab[:brk_pos+1]
                                daysc1 = days[:brk_pos+1]

                                mx = np.max(LAIc1)
                                mn = np.min(LAIc1)
                                thrs = mn + 0.2 * ( mx - mn )
                                if mx < 0.01:
                                        Onsets[0,nlat,nlon] = 0
                                        Offsets[0,nlat,nlon] = 0
                                else:
                                        tr = np.where(LAIc1 >= thrs)
                                        trar = tr[0]
                                        if len(trar) == len(LAIc1):
                                                Onsets[0,nlat,nlon] = 0
                                                Offsets[0,nlat,nlon] = 0
                                        else:
                                                kont = np.zeros(shape=(len(trar)))
                                                for i in range(len(trar)-1):
                                                        if trar[i+1]-trar[i] <= 2:
                                                                kont[i] = 1
                                                        elif trar[i+1]-trar[i] <= 2:
                                                                kont[i] = 1
                                                        else:
                                                                kont[i] = 0
                                                if trar[len(trar)-1]-trar[len(trar)-2] <= 2:
                                                        kont[len(trar)-1] = 1
                                                else:
                                                        kont[len(trar)-1] = 0
						
						kont2 = 0
                                                for i in range(len(kont)):
                                                        if kont[i] == 0:
                                                                kont2 = kont2 + 1
                                                                break
                                                        else:
                                                                kont2 = kont2 + 1

                                                if kont2 == len(trar):
                                                        trar2 = np.zeros(shape=(kont2))
                                                        for i in range(kont2):
                                                                trar2[i] = trar[i]
                                                elif kont[0] == 0:
                                                        trar2 = np.zeros(shape=(len(kont)-1))
                                                        for i in range(len(kont)-1):
                                                                trar2[i] = trar[i+1]
                                                else:
                                                        valbc2=0
                                                        valac2=0
                                                        for i in range(len(kont)):
                                                                if i <= kont2-1:
                                                                        valbc2 = valbc2+1
                                                                else:
                                                                        valac2 = valac2+1
                                                        if valbc2 > valac2:
                                                                trar2 = np.zeros(shape=(kont2))
                                                                for i in range(len(trar2)):
                                                                        trar2[i] = trar[i]
                                                        else:
                                                                trar2 = np.zeros(shape=(valac2))
                                                                for i in range(len(trar2)):
                                                                        trar2[i] = trar[i+kont2]

						Onset = trar2[0]
                                                Offset = trar2[len(trar2)-1]

                                                Onsets[0,nlat,nlon] = daysc1[int(Onset)]
                                                Offsets[0,nlat,nlon] = daysc1[int(Offset)]

			# analyse the second cycle as a shorter type 1 cycle
                                LAIc2 = Ab[brk_pos-1:]
                                daysc2 = days[brk_pos-1:]

                                mx = np.max(LAIc2)
                                mn = np.min(LAIc2)
                                thrs = mn + 0.2 * ( mx - mn )
                                if mx < 0.01:
                                        Onsets[1,nlat,nlon] = 0
                                        Offsets[1,nlat,nlon] = 0
                                else:
                                        tr = np.where(LAIc2 >= thrs)
                                        trar = tr[0]
                                        if len(trar) == len(LAIc2):
                                                Onsets[1,nlat,nlon] = 0
                                                Offsets[1,nlat,nlon] = 0
					else:
                                                kont = np.zeros(shape=(len(trar)))
                                                for i in range(len(trar)-1):
                                                        if trar[i+1]-trar[i] <= 2:
                                                                kont[i] = 1
                                                        elif trar[i+1]-trar[i] <= 2:
                                                                kont[i] = 1
                                                        else:
                                                                kont[i] = 0
                                                if trar[len(trar)-1]-trar[len(trar)-2] <= 2:
                                                        kont[len(trar)-1] = 1
                                                else:
                                                        kont[len(trar)-1] = 0

                                                kont2 = 0
                                                for i in range(len(kont)):
                                                        if kont[i] == 0:
                                                                kont2 = kont2 + 1
                                                                break
                                                        else:
                                                                kont2 = kont2 + 1

                                                if kont2 == len(trar):
                                                        trar2 = np.zeros(shape=(kont2))
                                                        for i in range(kont2):
                                                                trar2[i] = trar[i]
                                                elif kont[0] == 0:
                                                        trar2 = np.zeros(shape=(len(kont)-1))
                                                        for i in range(len(kont)-1):
                                                                trar2[i] = trar[i+1]
						else:
                                                        valbc2=0
                                                        valac2=0
                                                        for i in range(len(kont)):
                                                                if i <= kont2-1:
                                                                        valbc2 = valbc2+1
                                                                else:
                                                                        valac2 = valac2+1
                                                        if valbc2 > valac2:
                                                                trar2 = np.zeros(shape=(kont2))
                                                                for i in range(len(trar2)):
                                                                        trar2[i] = trar[i]
                                                        else:
                                                                trar2 = np.zeros(shape=(valac2))
                                                                for i in range(len(trar2)):
                                                                        trar2[i] = trar[i+kont2]

                                                Onset = trar2[0]
                                                Offset = trar2[len(trar2)-1]

                                                Onsets[1,nlat,nlon] = daysc2[int(Onset)]
                                                Offsets[1,nlat,nlon] = daysc2[int(Offset)]
		# check on the order of cycles, in case invert them
                                if not np.isnan(Onsets[0,nlat,nlon]) and not np.isnan(Onsets[1,nlat,nlon]):
                                        if Onsets[0,nlat,nlon] > Onsets[1,nlat,nlon]:
                                                tmpN0 = Onsets[0,nlat,nlon]
                                                tmpN1 = Onsets[1,nlat,nlon]
                                                tmpF0 = Offsets[0,nlat,nlon]
                                                tmpF1 = Offsets[1,nlat,nlon]
                                                Onsets[0,nlat,nlon] = tmpN1
                                                Onsets[1,nlat,nlon] = tmpN0
                                                Offsets[0,nlat,nlon] = tmpF1
                                                Offsets[1,nlat,nlon] = tmpF0

                        elif SeasType[nlat,nlon] == 3 and SeasClss[nlat,nlon] == 4:

			# determine the two cycles: one in the central 13 points and one using the 13 lateral points
                                LAIc1 = Ab[6:19]
                                daysc1 = days[6:19]

                                LAIc2_a = Ab[:7]
                                LAIc2_b = Ab[18:]
				daysc2_a = days[:7]
				daysc2_b = days[18:]
                                LAIc2 = np.zeros(shape=(len(LAIc2_a)+len(LAIc2_b)))
                                daysc2 = np.zeros(shape=(len(LAIc2_a)+len(LAIc2_b)))

				for ii in range(len(LAIc2_b)):
				        LAIc2[ii] = LAIc2_b[ii]
				        daysc2[ii] = daysc2_b[ii]
				
				for ii in range(len(LAIc2)-len(LAIc2_b)):
				        jj = len(LAIc2_b) + ii
				        LAIc2[jj] = LAIc2_a[ii]
				        daysc2[jj] = daysc2_a[ii]

#                                jj = 0
#                                for ii in range(len(Ab)):
#                                        if ii < 7 or ii >= 18:
#                                                LAIc2[jj] = Ab[ii]
#                                                daysc2[jj] = days[ii]
#                                                jj = jj+1

			# analyse the first cycle as a shorter type 1 cycle
                                mx = np.max(LAIc1)
                                mn = np.min(LAIc1)
                                thrs = mn + 0.2 * ( mx - mn )
                                if mx < 0.01:
                                        Onsets[0,nlat,nlon] = 0
                                        Offsets[0,nlat,nlon] = 0
                                else:
                                        tr = np.where(LAIc1 >= thrs)
                                        trar = tr[0]
                                        if len(trar) == len(LAIc1):
                                                Onsets[0,nlat,nlon] = 0
                                                Offsets[0,nlat,nlon] = 0
                                        else:
                                                kont = np.zeros(shape=(len(trar)))
                                                for i in range(len(trar)-1):
                                                        if trar[i+1]-trar[i] <= 2:
                                                                kont[i] = 1
                                                        elif trar[i+1]-trar[i] <= 2:
                                                                kont[i] = 1
                                                        else:
                                                                kont[i] = 0
                                                if trar[len(trar)-1]-trar[len(trar)-2] <= 2:
                                                        kont[len(trar)-1] = 1
                                                else:
                                                        kont[len(trar)-1] = 0
						
						kont2 = 0
                                                for i in range(len(kont)):
                                                        if kont[i] == 0:
                                                                kont2 = kont2 + 1
                                                                break
                                                        else:
                                                                kont2 = kont2 + 1

                                                if kont2 == len(trar):
                                                        trar2 = np.zeros(shape=(kont2))
                                                        for i in range(kont2):
                                                                trar2[i] = trar[i]
                                                elif kont[0] == 0:
                                                        trar2 = np.zeros(shape=(len(kont)-1))
                                                        for i in range(len(kont)-1):
                                                                trar2[i] = trar[i+1]
                                                else:
                                                        valbc2=0
                                                        valac2=0
                                                        for i in range(len(kont)):
                                                                if i <= kont2-1:
                                                                        valbc2 = valbc2+1
                                                                else:
                                                                        valac2 = valac2+1
                                                        if valbc2 > valac2:
                                                                trar2 = np.zeros(shape=(kont2))
                                                                for i in range(len(trar2)):
                                                                        trar2[i] = trar[i]
                                                        else:
                                                                trar2 = np.zeros(shape=(valac2))
                                                                for i in range(len(trar2)):
                                                                        trar2[i] = trar[i+kont2]

                                                Onset = trar2[0]
                                                Offset = trar2[len(trar2)-1]
				
						Onsets[0,nlat,nlon] = daysc1[int(Onset)]
                                                Offsets[0,nlat,nlon] = daysc1[int(Offset)]

                        # analyse the second cycle as a shorter type 1 cycle
				mx = np.max(LAIc2)
                                mn = np.min(LAIc2)
                                thrs = mn + 0.2 * ( mx - mn )
                                if mx < 0.01:
                                        Onsets[1,nlat,nlon] = 0
                                        Offsets[1,nlat,nlon] = 0
                                else:
                                        tr = np.where(LAIc2 >= thrs)
                                        trar = tr[0]
                                        if len(trar) == len(LAIc2):
                                                Onsets[1,nlat,nlon] = 0
                                                Offsets[1,nlat,nlon] = 0
                                        else:
                                                kont = np.zeros(shape=(len(trar)))
                                                for i in range(len(trar)-1):
                                                        if trar[i+1]-trar[i] <= 2:
                                                                kont[i] = 1
                                                        elif trar[i+1]-trar[i] <= 2:
                                                                kont[i] = 1
                                                        else:
                                                                kont[i] = 0
                                                if trar[len(trar)-1]-trar[len(trar)-2] <= 2:
                                                        kont[len(trar)-1] = 1
                                                else:
                                                        kont[len(trar)-1] = 0
						
						kont2 = 0
                                                for i in range(len(kont)):
                                                        if kont[i] == 0:
                                                                kont2 = kont2 + 1
                                                                break
                                                        else:
                                                                kont2 = kont2 + 1

                                                if kont2 == len(trar):
                                                        trar2 = np.zeros(shape=(kont2))
                                                        for i in range(kont2):
                                                                trar2[i] = trar[i]
                                                elif kont[0] == 0:
                                                        trar2 = np.zeros(shape=(len(kont)-1))
                                                        for i in range(len(kont)-1):
                                                                trar2[i] = trar[i+1]
                                                else:
                                                        valbc2=0
                                                        valac2=0
                                                        for i in range(len(kont)):
                                                                if i <= kont2-1:
                                                                        valbc2 = valbc2+1
                                                                else:
                                                                        valac2 = valac2+1
                                                        if valbc2 > valac2:
                                                                trar2 = np.zeros(shape=(kont2))
                                                                for i in range(len(trar2)):
                                                                        trar2[i] = trar[i]
                                                        else:
                                                                trar2 = np.zeros(shape=(valac2))
                                                                for i in range(len(trar2)):
                                                                        trar2[i] = trar[i+kont2]

                                                Onset = trar2[0]
                                                Offset = trar2[len(trar2)-1]

						Onsets[1,nlat,nlon] = daysc2[int(Onset)]
                                                Offsets[1,nlat,nlon] = daysc2[int(Offset)]

                # check on the order of cycles, in case invert them
                                if not np.isnan(Onsets[0,nlat,nlon]) and not np.isnan(Onsets[1,nlat,nlon]):
                                        if Onsets[0,nlat,nlon] > Onsets[1,nlat,nlon]:
                                                tmpN0 = Onsets[0,nlat,nlon]
                                                tmpN1 = Onsets[1,nlat,nlon]
                                                tmpF0 = Offsets[0,nlat,nlon]
                                                tmpF1 = Offsets[1,nlat,nlon]
                                                Onsets[0,nlat,nlon] = tmpN1
                                                Onsets[1,nlat,nlon] = tmpN0
                                                Offsets[0,nlat,nlon] = tmpF1
                                                Offsets[1,nlat,nlon] = tmpF0

		else:
			SeasType[nlat,nlon] = np.nan
			Onsets[0,nlat,nlon] = np.nan
			Onsets[1,nlat,nlon] = np.nan
			Offsets[0,nlat,nlon] = np.nan
			Offsets[1,nlat,nlon] = np.nan

# ******** Print netcdf output files ******

# **** create file
out_nc1 = Dataset(outname1, "w", format="NETCDF4")

# **** define dimensions
lat = out_nc1.createDimension('lat', len(lats))
lon = out_nc1.createDimension('lon', len(lons))
pks = out_nc1.createDimension('peak', 2) 

# **** define varaibles
latitudes = out_nc1.createVariable('lat', np.float32, ('lat',))
longitudes = out_nc1.createVariable('lon', np.float32, ('lon',))

Seas_Ampl = out_nc1.createVariable('Seas_Ampl',np.float32,('lat','lon'))
onset = out_nc1.createVariable('onset',np.float32,('peak','lat','lon'))
offset = out_nc1.createVariable('offset',np.float32,('peak','lat','lon'))
Seas_Type = out_nc1.createVariable('Seas_Type',np.float32,('lat','lon'))
Seas_Class = out_nc1.createVariable('Seas_Class',np.float32,('lat','lon'))
Seas_Class1 = out_nc1.createVariable('Seas_Class1',np.float32,('lat','lon'))

# **** define global attribute
out_nc1.description = 'Contain growing season onset, offset and seasonal amplitude '
out_nc1.source = 'LAI3g'
out_nc1.history = 'Created ' + str(datetime.now())

# **** define variables attribute
latitudes.long_name = 'latitude'
latitudes.units = 'degree_north'

longitudes.long_name = 'longitude'
longitudes.units = 'degree_east'

Seas_Ampl.long_name = 'LAI Seasonal Amplitude'
Seas_Ampl.units = 'm2/m2'

onset.long_name = 'Growing Season Onset'
onset.units = 'm2/m2'

offset.long_name = 'Growing Season Offset'
offset.units = 'm2/m2'

Seas_Type.long_name = 'LAI Seasonal Type'
Seas_Type.units = 'index'

Seas_Class.long_name = 'Two peaks Seasonal Class'
Seas_Class.units = 'index'

Seas_Class1.long_name = 'One peak Seasonal Class'
Seas_Class1.units = 'index'

# **** assign values to variables
latitudes[:] = lats
longitudes[:] = lons

Seas_Ampl[:,:] = SeasAmpl 
onset[:,:,:] = Onsets
offset[:,:,:] = Offsets
Seas_Type[:,:] = SeasType 
Seas_Class[:,:] = SeasClss
Seas_Class1[:,:] = SeasClss1

out_nc1.close()

print 'Writing:  ', outname1

