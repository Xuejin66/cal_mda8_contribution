import numpy as np
import netCDF4 as nc

####################################更新备注#################################
# 1.16 完成数据提取模块，完成维数转化模块，完成MDA8计算模块,完成MDA8索引提取模块

def prep_var(path, filename):
    data = nc.Dataset(path + filename)
    varnames = list(data.variables.keys())  #read the nc variable name as a list
    var = data['O3_R1R2R3'][23:,0,:,:]
    data.close()
    return var

def trans2var4d(var):
    var_4d = np.zeros((31,24,216,192), dtype=float)
    for day_4d in range(31):
        for hour_4d in range(24):
            var_4d[day_4d,hour_4d,:,:] = var[24*day_4d + hour_4d,:,:]
    return var_4d

def O3_4d2nc(var,name):
    fout = nc.Dataset(name,'w',format = 'NETCDF4')
    fout.createDimension('day',31)  
    fout.createDimension('hour',24) 
    fout.createDimension('row',216)   
    fout.createDimension('col',192)
    fout.createVariable('day',np.int,('day'))  
    fout.createVariable('hour',np.int,('hour'))  
    fout.createVariable('row',np.int,('row'))  
    fout.createVariable('col',np.int,('col'))
    day = np.arange(31)
    hour = np.arange(24)
    row = np.arange(216)
    col = np.arange(192)
    fout.variables['day'][:] = day
    fout.variables['hour'][:] = hour
    fout.variables['row'][:] = row
    fout.variables['col'][:] = col
    fout.createVariable('O3_4d', np.float32, ('day','hour','row','col'))
    fout.variables['O3_4d'][:] = var
    fout.close()

def cal_mda8(var_4d):
    mda8 = np.zeros((31,216,192),dtype = float)
    mda8_TimeInd = np.zeros((31,24,216,192),dtype = int)
    for col in range(192):
        for row in range(216):
            for day in range(31):
                #print('NO.%d day'%day)
                mda_8 = 0
                mda8_list = []
                for hour in range(16):
                    #print('NO.%d hour'%hour)
                    d8 = var_4d[day, hour:hour + 8, row, col]
                    da8 = np.nanmean(d8)
                    mda8_list.append(da8)                
                    #if da8 >= mda_8:
                        #mda_8 = da8
                mda8[day,row,col] = max(mda8_list)
                time_index = mda8_list.index(max(mda8_list))
                mda8_TimeInd[day,time_index:time_index + 8,row,col] = 1
    return mda8,mda8_TimeInd

def mda8_2nc(var,name):
    fout = nc.Dataset(name,'w',format = 'NETCDF4')
    fout.createDimension('day',31)  
    fout.createDimension('row',216)   
    fout.createDimension('col',192)
    fout.createVariable('day',np.int,('day'))   
    fout.createVariable('row',np.int,('row'))  
    fout.createVariable('col',np.int,('col'))
    day = np.arange(31)
    row = np.arange(216)
    col = np.arange(192)
    fout.variables['day'][:] = day
    fout.variables['row'][:] = row
    fout.variables['col'][:] = col
    fout.createVariable( 'mda8', np.float32, ('day','row','col'))
    fout.variables['mda8'][:] = var
    fout.close()

def mda8TimeInd2nc(var,name):
    fout = nc.Dataset(name,'w',format = 'NETCDF4')
    fout.createDimension('day',31)  
    fout.createDimension('hour',24) 
    fout.createDimension('row',216)   
    fout.createDimension('col',192)
    fout.createVariable('day',np.int,('day'))  
    fout.createVariable('hour',np.int,('hour'))  
    fout.createVariable('row',np.int,('row'))  
    fout.createVariable('col',np.int,('col'))
    day = np.arange(31)
    hour = np.arange(24)
    row = np.arange(216)
    col = np.arange(192)
    fout.variables['day'][:] = day
    fout.variables['hour'][:] = hour
    fout.variables['row'][:] = row
    fout.variables['col'][:] = col
    fout.createVariable('mda8TimeInd', np.int32, ('day','hour','row','col'))
    fout.variables['mda8TimeInd'][:] = var
    fout.close()

path = '/data3/liuhanqing/projects/postproc/2.combine/output/201807_OSAT/'
filename = 'camx.YRD4km.R1R2R3.201807.OSAT.sa.grd01.ncf'

O3_3d = prep_var(path,filename)
O3_4d = trans2var4d(O3_3d)
mda8,mda8_TimeInd = cal_mda8(O3_4d)
print(mda8_TimeInd[6,:,125,125])
O3_4d2nc(O3_4d,'O3_4d.nc')
mda8_2nc(mda8,'mda8.nc')
mda8TimeInd2nc(mda8_TimeInd,'mda8_TimeInd.nc')

###########################################
#         DATA VERIFICATION               #
###########################################
print('July 7th O3_hour')
print(O3_4d[6,:,125,125])
print('MDA8 In July')
print(mda8[:,125,125])
print('The time index of MDA8 in July 7th')

