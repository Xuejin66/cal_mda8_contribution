#!/usr/people/zliu/miniconda/bin/python

from scipy.io import netcdf
import numpy as np
# from mpl_toolkits.basemap import Basemap
# from netCDF4 import Dataset
import math
import sys


# read netcdf file with using the filename and the variable name
def get_netcdf3_var(ncdf_file, varname):
    data = netcdf.netcdf_file(ncdf_file, 'r')
    var = data.variables[varname]
    return var

def get_tracernames(spec):
    #ext = ['BC','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O',\
    ext = ['BC','A','B','C','D','E','BY','FY','IC']
    ntrac = len(ext)
    tracernames = [spec+'_'+ext[i] for i in range(ntrac)]
    tracernames[0]=spec+ext[0];tracernames[-1]=spec+ext[-1] # IC, BC
    return tracernames


def read_var3d_hrly ( ncdf_file , spec ) :
    var3d = get_netcdf3_var ( ncdf_file , spec ).data[:,0,:,:] 
    return var3d

def get_latlon ( grid_file ) :
    lats_ctr = get_netcdf3_var ( grid_file , 'LAT' )[0,0,:,:]
    lons_ctr = get_netcdf3_var ( grid_file , 'LON' )[0,0,:,:]
    return lats_ctr , lons_ctr


def get_dimensions ( var ) :
    if var.ndim == 4 :
       nx = var.shape[-1]; ny = var.shape[-2]; nsteps = var.shape[-4]
    if var.ndim == 3 :
       nx = var.shape[-1]; ny = var.shape[-2]; nsteps = var.shape[-3]
    return nx , ny , nsteps

def calcvar_8hra_1day ( var3d_hrly_48hr ) :
    # returns running averages for 24 hours
    # input : hourly averages for 48 hours [48,ny,nx]; 
    #         48 hours is enough to allow time shifting. 
    # output: 8-hour averages for 24 hours [24,ny,nx]
    nx , ny , nsteps = get_dimensions ( var3d_hrly_48hr )
    var3d_8hra_1day = np.zeros ( (24 , ny , nx) )
    for t in range ( 24 ):
        var3d_8hra_1day[t,:,:] = np.average(var3d_hrly_48hr[t:t+8,:,:],axis=0) # this cannot be k:k+7
    return var3d_8hra_1day

def calcvar_1hrmax_1day ( var3d_hrly_24hr ) :
    # returns hourly maximum for 24 hours
    # input : hourly conc for 24 hours [24,ny,nx];
    # output: maximum of 24 hours [ny,nx]
    nx , ny , nsteps = get_dimensions ( var3d_hrly_24hr )
    var2d_1hrmax_1day = np.amax( var3d_hrly_24hr, axis = 0)
    return var2d_1hrmax_1day

def calcvar_mda8_ind_1day ( var3d_8hra_24hr ) :
    # returns MDA8 and its hour index
    # input : 8-hour running averages for 24 hours  [24,ny,nx]
    # output: MDA8  [ny,nx]; index [ny,nx]
    var2d_mda8_1day = np.amax ( var3d_8hra_24hr , axis = 0 )
    ind2d_mda8_1day = np.argmax ( var3d_8hra_24hr , axis=0 )
    return var2d_mda8_1day , ind2d_mda8_1day

def getvar_mda8hr_1day ( var3d_8hra_24hr , ind2d_mda8_1day ) :
    # returns (tracer) values at the hour when MDA8 occurs
    # input : 8-hour running averages for 24 hours  [24,ny,nx]
    #         MDA8-hour index [24,ny,nx] 
    nx , ny , nsteps = get_dimensions ( var3d_8hra_24hr ) 
    var2d_mda8hr_1day = np.zeros((ny,nx))
    for j in range ( ny ) :
        for i in range ( nx ) :
            var2d_mda8hr_1day[j,i]=var3d_8hra_24hr[ind2d_mda8_1day[j,i],j,i] 
    return var2d_mda8hr_1day

def write_output_netcdf_for_ncl ( spec , var , outfile) :
    from scipy.io.netcdf import NetCDFFile as Dataset
    from numpy import arange, dtype # array module from http://numpy.scipy.org

    if var.ndim == 3: # it is a 3d variable
	        nx = var.shape[-1]; ny = var.shape[-2]; nstep = var.shape[-3];
	        ncfile = Dataset(outfile,'w') 
		# create the dimensions.
	        ncfile.createDimension('nx',nx)
	        ncfile.createDimension('ny',ny)
	        ncfile.createDimension('nstep',nstep)
	        varout = ncfile.createVariable(spec,dtype('float32').char,('nstep','ny','nx'))
	        varout.units = 'ug/m3'
	        varout[:,:,:] = var[:,:,:]
		# close the file.
	        ncfile.close()
    if var.ndim == 2: # it is a 2d variable
        nx = var.shape[-1]; ny = var.shape[-2];
        ncfile = Dataset(outfile,'w')
        # create the dimensions.
        ncfile.createDimension('nx',nx)
        ncfile.createDimension('ny',ny)
        varout = ncfile.createVariable(spec,dtype('float32').char,('ny','nx'))
        varout.units = 'ug/m3'
        varout[:,:] = var[:,:]
        # close the file.
        ncfile.close() 
        print ('*** SUCCESS writing netCDF file')
    return 

def calc_mda8 ( ncf_file, spec , tshift ):

    # read the hourly concentrations 
    var3d_hrly = read_var3d_hrly ( ncf_file , spec  )
    
    # get the dimensions
    nx , ny , nsteps = get_dimensions ( var3d_hrly )
    ndays =nsteps // 24
    # define the variables
    var3d_mda8_ndays = np.zeros((ndays,ny,nx))
    ind3d_mda8_ndays = np.zeros((ndays,ny,nx))  
    
    for d in range ( ndays ) :
        step0 = d * 24 + tshift         
        var3d_hrly_48hr = var3d_hrly[step0:step0+48,:,:]
        var3d_8hra_24hr = calcvar_8hra_1day ( var3d_hrly_48hr )
        var2d_mda8_1day , ind2d_mda8_1day = calcvar_mda8_ind_1day ( var3d_8hra_24hr ) 
        print ('calculate MDA8 O3 for day : ' , d+1 , '(10,10) = ' , var2d_mda8_1day[10,10], 'hour = ', ind2d_mda8_1day[10,10])
        var3d_mda8_ndays[d,:,:] = var2d_mda8_1day[:,:]
        ind3d_mda8_ndays[d,:,:] = ind2d_mda8_1day[:,:]
    return var3d_mda8_ndays , ind3d_mda8_ndays         

def calcvar_h4mda8 ( var3d_mda8_ndays ) :
    nx , ny , nsteps = get_dimensions ( var3d_mda8_ndays )
    ind3d_mda8_ndays_sorted = np.argsort ( var3d_mda8_ndays , axis = 0 )
    ind2d_h4mda8 = ind3d_mda8_ndays_sorted[-4,:,:]
    var2d_h4mda8 = np.zeros((ny,nx))
    for i in range(nx) :
        for j in range(ny) :
            var2d_h4mda8[j,i] = var3d_mda8_ndays[ind2d_h4mda8[j,i],j,i]    
    return var2d_h4mda8 , ind2d_h4mda8

def calcvar_avg_mda8 ( var3d_mda8_ndays ) :
    nx , ny , nsteps = get_dimensions ( var3d_mda8_ndays )
    var2d_avg_mda8 = np.average(var3d_mda8_ndays,axis = 0)
    print (var2d_avg_mda8[10,10])
    return var2d_avg_mda8

def calcvar_annavg_1hrmax ( ncf_file, spec , tshift ):
        var3d_hrly_24hr = var3d_hrly[step0:step0+24,:,:]
        # read the hourly concentrations
        var3d_hrly = read_var3d_hrly ( ncf_file , spec  )
  
        # get the dimensions
        nx , ny , nsteps = get_dimensions ( var3d_hrly )
        ndays = nsteps / 24 - 1
        # define the variables
        var3d_1hrmax_ndays = np.zeros((ndays,ny,nx))
        ind3d_1hrmax_ndays = np.zeros((ndays,ny,nx))

        for d in range ( ndays ) :
            step0 = d * 24 + tshift
            var3d_hrly_24hr = var3d_hrly[step0:step0+24,:,:]
            var2d_1hrmax_1day = calcvar_1hrmax_1day( var3d_hrly_24hr )
            var3d_1hrmax_ndays[d,:,:] = var2d_1hrmax_1day[:,:]
        var2d_annavg_1hrmax = np.average(var3d_1hrmax_ndays, axis = 0 )      
        return var2d_annavg_1hrmax

def calcvar_90thmda8 ( var3d_mda8_ndays ) :
    nx , ny , nsteps = get_dimensions ( var3d_mda8_ndays )
    ind3d_mda8_ndays_sorted = np.argsort ( var3d_mda8_ndays , axis = 0 )
    ind = int(nsteps*0.1)
    print (ind)
    ind2d_90thmda8 = ind3d_mda8_ndays_sorted[-ind,:,:] # 90th highest mda8
    var2d_90thmda8 = np.zeros((ny,nx))
    for i in range(nx) :
        for j in range(ny) :
            var2d_90thmda8[j,i] = var3d_mda8_ndays[ind2d_90thmda8[j,i],j,i]
    return var2d_90thmda8 , ind2d_90thmda8


def calcvar_dailyavg_ndays ( ncdf_file , spec, tshift ) :
	var3d_hrly = read_var3d_hrly ( ncf_file , spec  )
	nx , ny , nsteps = get_dimensions ( var3d_hrly )
	ndays = nsteps / 24 
	var3d_dailyavg_ndays = np.zeros((ndays,ny,nx))
	for d in range ( ndays ) :
		step0 = d * 24 + tshift
		var3d_hrly_24hr = var3d_hrly[step0:step0+24,:,:]
		var3d_dailyavg_ndays[d,:,:] = np.average(var3d_hrly_24hr, axis=0)
	return var3d_dailyavg_ndays

def calcvar_1hrmax_ndays ( ncdf_file , spec, tshift ) :
    var3d_hrly = read_var3d_hrly ( ncf_file , spec  )
    nx , ny , nsteps = get_dimensions ( var3d_hrly )
    ndays = nsteps / 24
    var3d_1hrmax_ndays = np.zeros((ndays,ny,nx))
    for d in range ( ndays ) :
        step0 = d * 24 + tshift
        var3d_hrly_24hr = var3d_hrly[step0:step0+24,:,:]
        var3d_1hrmax_ndays[d,:,:] = np.max(var3d_hrly_24hr, axis=0)
    return var3d_1hrmax_ndays


def calcvar_daytimeavg_ndays ( ncdf_file , spec, tshift, hr0, hr1 ) :
    var3d_hrly = read_var3d_hrly ( ncf_file , spec  )
    nx , ny , nsteps = get_dimensions ( var3d_hrly )
    ndays = nsteps / 24
    var3d_daytimeavg_ndays = np.zeros((ndays,ny,nx))
    for d in range ( ndays ) :
        step0 = d * 24 + tshift + hr0
        step1 = d * 24 + tshift + hr1 + 1
        var3d_daytimeavg_ndays[d,:,:] = np.average(var3d_hrly[step0:step1,:,:], axis=0)
    return var3d_daytimeavg_ndays

def calcvar_epiavg ( ncdf_file, spec, tshift ):
    var3d_hrly = read_var3d_hrly ( ncf_file , spec  )
    var2d_epiavg = np.average(var3d_hrly, axis=0)
    return var2d_epiavg	

def calcvar_epimax ( ncdf_file, spec, tshift ):
    var3d_hrly = read_var3d_hrly ( ncf_file , spec  )
    var2d_epimax = np.max(var3d_hrly, axis=0)
    return var2d_epimax

#####################################################################################
########## ONLY NEED TO MODIFY THIS SECTION ##############################

# !!!!! input files are always hourly data generated by combine step in Beijing time
#dir = '/data3/wangq/projects/postproc/2.combine/output/YRD4km/2018/07/MEGAN31.CMAQICBC/'
dir = '/data3/liuhanqing/projects/postproc/2.combine/output/201807_OSAT/'
ncf_file = dir+'camx.YRD4km.201807.OSAT.sa.grd01.ncf.BST'

tshift  = 0      # always prepare input data in Beijing time so that tshift = 0
spec    = 'O3_R1G9'
case    = 'Suzhou_OSAT_O3_R1G9'   # set label for your case
hr0 	= 10 	 # set starting hour if need to calculate average for selected times
hr1 	= 14 	 # set ending hour if need to calculate average for selected times
metric  = 'mda8_byday' 	
outdir = dir
						# choose from 'mda8_byday' daily mda8 for ndays
                 		# 			  'mda8_avg' averaged mda8
               		  	# 			  '90thmda8' 90th mda8
                 		# 			  'dailyavg_byday' daily average values for ndays
                 		# 	 		  '1hrmax_byday' daily max values for ndays
                 		# 	 		  '1hrmax' averaged daily max values
                 		# 			  'daytimeavg_byday' day time average depending on user defined day time range
                 		# 			  'daytimeavg' averages for selected day time average
  				 		# 			  'epiavg'   episode average (whether it is monthly average or annual average)
                 		# 			  'epimax'   episode maximum

# define output file
outfile = outdir + spec+"."+case+"."+metric+".nc"

##########################################################################################
########### DO NOT MODIFY BELOW THIS LINE UNLESS YOU KNOW WHAT TO DO ########################

if metric == 'mda8_byday':
	if spec != 'O3_R1G9':
		print ('***** ERROR ******')
		print (metric+' only works for O3_R1G9')
	else:
		print ('Calculating daily MDA8 O3_R1G9 concentrations')
		var3d_mda8_ndays , ind3d_mda8_ndays = calc_mda8 ( ncf_file, spec, tshift )
		write_output_netcdf_for_ncl ( spec , var3d_mda8_ndays , outfile)     

if metric == 'mda8_avg':
	if spec != 'O3':
		print ('**** ERROR ******')
		print (metric+' only works for O3')
	else:
		print ('Calculating episode average MDA8 O3 concentrations')
		var3d_mda8_ndays , ind3d_mda8_ndays = calc_mda8 ( ncf_file, spec, tshift )
		var2d_avg_mda8 = calcvar_avg_mda8(var3d_mda8_ndays)
		write_output_netcdf_for_ncl ( spec , var2d_avg_mda8 , outfile)     

if metric == '90thmda8':
    if spec != 'O3':
        print ('**** ERROR ******')
        print (metric+' only works for O3')
    else:
                print ('Calculating episode average MDA8 O3 concentrations')
                var3d_mda8_ndays , ind3d_mda8_ndays = calc_mda8  ( ncf_file, spec, tshift )
                var2d_90thmda8 , ind2d_90thmda8 = calcvar_90thmda8 ( var3d_mda8_ndays)
                write_output_netcdf_for_ncl ( spec , var2d_90thmda8 , outfile)     

if metric == 'dailyavg_byday':
        print ('Calculating daily average values for ' + spec)
        var3d_dailyavg_ndays = calcvar_dailyavg_ndays (ncf_file, spec, tshift )
        write_output_netcdf_for_ncl ( spec , var3d_dailyavg_ndays , outfile)     

if metric == '1hrmax_byday':
      print ('Calculating daily maximum values for ' + spec)
      var3d_1hrmax_ndays = calcvar_1hrmax_ndays (ncf_file, spec, tshift )
      write_output_netcdf_for_ncl ( spec , var3d_1hrmax_ndays , outfile)

if metric == '1hrmax':
    print ('Calculating averaged daily maximum values for ' + spec)

    var3d_1hrmax_ndays = calcvar_1hrmax_ndays (ncf_file, spec, tshift )
    var2d_1hrmax = np.average(var3d_1hrmax_ndays, axis=0)
    write_output_netcdf_for_ncl ( spec , var2d_1hrmax , outfile)

if metric == 'daytimeavg_byday':
        print ('Calculating average values for ' + str(hr0) + ':00 to ' + str(hr1) + ':00 for ' + spec)
        
        if hr0 > hr1:
                print ('****ERROR**** starting hour could not be greater than ending hour')
        else:
                var3d_daytimeavg_ndays = calcvar_daytimeavg_ndays ( ncf_file, spec, tshift, hr0, hr1 )
        write_output_netcdf_for_ncl ( spec , var3d_daytimeavg_ndays , outfile)

if metric == 'daytimeavg':
        print ('Calculating average values for ' + str(hr0) + ':00 to ' + str(hr1) + ':00 for ' + spec)
        if hr0 > hr1:
                print ('****ERROR**** starting hour could not be greater than ending hour')
        else:
                var3d_daytimeavg_ndays = calcvar_daytimeavg_ndays ( ncf_file, spec, tshift, hr0, hr1 )
                var2d_daytimeavg = np.average( var3d_daytimeavg_ndays , axis = 0 ) 
        write_output_netcdf_for_ncl ( spec , var2d_daytimeavg , outfile)

if metric == 'epiavg':
        print ('Calculating episode average values for ' + spec)
        var2d_epiavg = calcvar_epiavg ( ncf_file, spec, tshift )
        write_output_netcdf_for_ncl ( spec , var2d_epiavg , outfile )

if metric == 'epimax':
        print ('Calculating episode maximum values for ' + spec)
        var2d_epimax = calcvar_epimax ( ncf_file, spec, tshift )
        write_output_netcdf_for_ncl ( spec , var2d_epimax , outfile )
