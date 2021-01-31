import numpy as np
import netCDF4 as nc
import datetime

##########################  FUNCTION EXPLAINATION  ###################################
#
#      treat_unSAfile()  提取总臭氧小时浓度数据，并计算MDA8和MDA8每日的时间索引         
#      mda8_2nc()        将计算得到的MDA8时间索引数据写入netCDF文件中                  
#      mda8Ind_2nc()     将计算得到的MDA8数据写入netCDF文件中              
#            
######################################################################################

###########################################
#           DEFINE FUNCTION               #
###########################################

def treat_unSAfile(unSAfilePath,unSAfilename,unSAvarname,nhour,nrow,ncol):
    unSAdata = nc.Dataset(unSAfilePath + unSAfilename)
    unSApol = unSAdata[unSAvarname][:nhour,0,:,:]
    unSAdata.close()
    nday = nhour // 24 - 1
    mda8list = np.zeros((nday,nrow,ncol),dtype = float)
    mda8Ind = np.zeros((nday,31,nrow,ncol),dtype = int)
    mda8Max = 0
    for row in range(nrow):
        for col in range(ncol):
            for hour in range(nhour - 24):
                day = hour // 24
                # print('第%s行,第%s列,第%s天'%(row,col,day))
                mda8 = np.nanmean(unSApol[hour:hour + 8,row,col])                
                if mda8 >= mda8Max:
                    mda8Max = mda8
                    mda8MaxInd = hour               
                if (hour + 1) % 24 == 0:
                    mda8Ind[day,mda8MaxInd%24:mda8MaxInd%24 + 8,row,col] = 1
                    mda8list[day,row,col] = mda8Max
                    mda8Max = 0               
    return mda8Ind,mda8list

def mda8_2nc(mda8list,nhour,nrow,ncol):
    fout = nc.Dataset('mda8.nc','w',format = 'NETCDF4')
    nday = nhour // 24 - 1
    fout.createDimension('day',nday)  
    fout.createDimension('row',nrow)   
    fout.createDimension('col',ncol)
    fout.createVariable('day',np.int,('day'))   
    fout.createVariable('row',np.int,('row'))  
    fout.createVariable('col',np.int,('col'))
    day = np.arange(nday)
    row = np.arange(nrow)
    col = np.arange(ncol)
    fout.variables['day'][:] = day
    fout.variables['row'][:] = row
    fout.variables['col'][:] = col
    fout.createVariable( 'mda8', np.float32, ('day','row','col'))
    fout.variables['mda8'][:] = mda8list
    fout.close()

def mda8Ind_2nc(mda8Ind,nhour,ncrow,ncol):
    fout = nc.Dataset('mda8Ind.nc','w',format = 'NETCDF4')
    nday = nhour // 24 - 1
    fout.createDimension('day',nday)  
    fout.createDimension('hour',31) 
    fout.createDimension('row',nrow)   
    fout.createDimension('col',ncol)
    fout.createVariable('day',np.int,('day'))  
    fout.createVariable('hour',np.int,('hour'))  
    fout.createVariable('row',np.int,('row'))  
    fout.createVariable('col',np.int,('col'))
    day = np.arange(nday)
    hour = np.arange(31)
    row = np.arange(nrow)
    col = np.arange(ncol)
    fout.variables['day'][:] = day
    fout.variables['hour'][:] = hour
    fout.variables['row'][:] = row
    fout.variables['col'][:] = col
    fout.createVariable('mda8Ind', np.int32, ('day','hour','row','col'))
    fout.variables['mda8Ind'][:] = mda8Ind
    fout.close()

###########################################
#              INPUT AREA                 #
###########################################

   # 这一行填写未源解析netCDF文件路径
unSAfilePath = '/data3/liuhanqing/projects/postproc/2.combine/output/201807_OSAT/'   
   # 这一行填写未源解析netCDF文件名
# unSAfilename = 'camx.YRD4km.R1R2R3_new.201807.OSAT.sa.grd01.ncf'
unSAfilename = 'camx.YRD4km.Suzhou_fin.201807.OSAT.sa.grd01.ncf'   
   # 总臭氧小时浓度的变量名                      
unSAvarname = 'O3_R1R2R3'  
   # 网格规格和源解析小时数
nhour = 744                                                          
nrow = 216
ncol = 192


###########################################
#              MAIN PROGRAM               #
###########################################

begin = datetime.datetime.now() 
print("START TIME: %s. It is estimated to take about 70 minutes, please be patient. " %begin)

  # 提取总臭氧小时浓度数据，并计算MDA8和MDA8每日的时间索引
mda8Ind,mda8list = treat_unSAfile(unSAfilePath,unSAfilename,unSAvarname,nhour,nrow,ncol)
  # 将计算得到的MDA8时间索引数据写入netCDF文件中
mda8Ind_2nc(mda8Ind,nhour,nrow,ncol)
  # 将计算得到的MDA8数据写入netCDF文件中
mda8_2nc(mda8list,nhour,nrow,ncol)

end = datetime.datetime.now()
print ("END TIME: %s. The program was completed successfully! " %end)