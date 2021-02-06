import numpy as np
import netCDF4 as nc
import datetime

############################################################  使用说明  #######################################################
#
#  你需要combine两个nc文件，一个nc文件中含有总臭氧浓度变量，即所有源贡献和,将该nc文件的路径和文件名输入在下方'O3TOT'相关地址中。
#  另一个文件为源解析输出文件，含有各类臭氧源贡献变量；将该nc文件的路径和文件名输入在下方'O3SA'相关地址中。
#
#################################################################################################################################

############################################################  待更新内容  #######################################################
#
#  1、优化标记变量模块，减少重复代码。
#  
#
#################################################################################################################################

def treat_O3TOTdata(O3TOTfilePath,O3TOTfilename,O3TOTvarname,nhour,nrow,ncol):
    O3TOTdata = nc.Dataset(O3TOTfilePath + O3TOTfilename)
    O3TOT = O3TOTdata[O3TOTvarname][:nhour,0,:,:]
    O3TOTdata.close()
    DA8TOT = np.zeros((nhour - 8,nrow,ncol),dtype = float)
    for hour in range(nhour - 8):
        DA8TOT[hour,...] = np.nanmean(O3TOT[hour:hour + 8],axis = 0)
    nday = nhour // 24 - 1
    DA8TOT_4d = np.zeros((nday,24,nrow,ncol), dtype=float)
    for day in range(nday):
        # for hour in range(31):
        DA8TOT_4d[day,...] = DA8TOT[24*day:24*(day + 1),...]  
    return DA8TOT_4d

def treat_O3SAdata_CAMx(O3SAfilePath,O3SAfileName,nsector,nhour,nrow,ncol):
    O3SAdata = nc.Dataset(O3SAfilePath + O3SAfileName)
    varnames = list(O3SAdata.variables.keys())  #read the nc variable name as a list
    varnames_new = []
    for varname in varnames:
        if 'O3' in varname[:2]:
            varnames_new.append(varname)
    O3SA_5d = np.zeros((nhour,2,nsector,nrow,ncol),dtype = float)  
    for i in range(2*nsector):
        O3SA_5d[:,i%2,i//2,:,:] = O3SAdata[varnames_new[i]] [:nhour,0,:,:]
    O3SAdata.close()   
    DA8SA = np.zeros((nhour - 8,2,nsector,nrow,ncol),dtype = float) 
    for hour in range(nhour - 8):
        DA8SA[hour,...] = np.nanmean(O3SA_5d[hour:hour + 8],axis = 0)
    nday = nhour // 24 - 1
    DA8SA_6d = np.zeros((nday,24,2,nsector,nrow,ncol), dtype=float)
    for day in range(nday):
        # for hour in range(31):
        DA8SA_6d[day,...] = DA8SA[24*day:24*(day + 1),...]
    return DA8SA_6d

def treat_O3SAdata_CMAQ(O3SAfilePath,O3SAfileName,nsector,nhour,nrow,ncol):
    O3SAdata = nc.Dataset(O3SAfilePath + O3SAfileName)
    varnames = list(O3SAdata.variables.keys())  #read the nc variable name as a list
    varnames_new = []
    for varname in varnames:
        if 'O3' in varname[:2]:
            varnames_new.append(varname)
    O3SA_5d = np.zeros((nhour,nsector,nrow,ncol),dtype = float)  
    for sector in range(nsector):
        O3SA_5d[:,sector,:,:] = O3SAdata[varnames_new[sector]] [:nhour,0,:,:]
    O3SAdata.close()   
    DA8SA = np.zeros((nhour - 8,nsector,nrow,ncol),dtype = float) 
    for hour in range(nhour - 8):
        DA8SA[hour,...] = np.nanmean(O3SA_5d[hour:hour + 8],axis = 0)
    nday = nhour // 24 - 1
    DA8SA_6d = np.zeros((nday,24,nsector,nrow,ncol), dtype=float)
    for day in range(nday):
        # for hour in range(31):
        DA8SA_6d[day,...] = DA8SA[24*day:24*(day + 1),...]
    return DA8SA_6d

def cal_MDA8SA_CAMx(DA8TOT,DA8SA_6d,nsector,nhour,nrow,ncol):
    MDA8_IND = np.argmax(DA8TOT,axis = 1)
    nday = nhour // 24 - 1
    MDA8SA = np.zeros((nday,2,nsector,nrow,ncol), dtype=float)
    for day in range(nday):
        for row in range(nrow):
            for col in range(ncol):
                for spec in range(2):
                    for sector in range(nsector):
                        TimeIndex = MDA8_IND[day,row,col]
                        MDA8SA[day,spec,sector,row,col] = DA8SA_6d[day,TimeIndex,spec,sector,row,col] 
    return MDA8SA

def cal_MDA8SA_CMAQ(DA8TOT,DA8SA_6d,nsector,nhour,nrow,ncol):
    MDA8_IND = np.argmax(DA8TOT,axis = 1)
    nday = nhour // 24 - 1
    MDA8SA = np.zeros((nday,nsector,nrow,ncol), dtype=float)
    for day in range(nday):
        for row in range(nrow):
            for col in range(ncol):
                for sector in range(nsector):
                    TimeIndex = MDA8_IND[day,row,col]
                    MDA8SA[day,sector,row,col] = DA8SA_6d[day,TimeIndex,sector,row,col] 
    return MDA8SA

def save_MDA8SA_CAMx(MDA8SA,Variables,VariableNames,nhour,nrow,ncol):
    nday = nhour // 24 - 1
    fout = nc.Dataset('MDA8SA.nc','w',format = 'NETCDF4')
    fout.createDimension('day',nday)  
    fout.createDimension('species',2)  
    fout.createDimension('sector',nsector)  
    fout.createDimension('row',nrow)   
    fout.createDimension('col',ncol)
    fout.createVariable('day',np.int,('day')) 
    fout.createVariable('species',np.int,('species'))
    fout.createVariable('sector',np.int,('sector'))  
    fout.createVariable('row',np.int,('row'))  
    fout.createVariable('col',np.int,('col'))
    day = np.arange(nday)
    species = np.arange(2)
    sector = np.arange(nsector)
    row = np.arange(nrow)
    col = np.arange(ncol)
    fout.variables['day'][:] = day
    fout.variables['species'][:] = species
    fout.variables['sector'][:] = sector
    fout.variables['row'][:] = row
    fout.variables['col'][:] = col
    fout.createVariable('MDA8SA', np.float32, ('day','species','sector','row','col'))
    fout.variables['MDA8SA'][:] = MDA8SA
    for num in range(len(Variables)):
        fout.createVariable(VariableNames[num], np.float32,('day','row','col'))
        fout.variables[VariableNames[num]][:] = Variables[num]
    fout.close()

def save_MDA8SA_CMAQ(MDA8SA,Variables,VariableNames,nhour,nrow,ncol):
    nday = nhour // 24 - 1
    fout = nc.Dataset('MDA8SA.nc','w',format = 'NETCDF4')
    fout.createDimension('day',nday)   
    fout.createDimension('sector',nsector)  
    fout.createDimension('row',nrow)   
    fout.createDimension('col',ncol)
    fout.createVariable('day',np.int,('day')) 
    fout.createVariable('sector',np.int,('sector'))  
    fout.createVariable('row',np.int,('row'))  
    fout.createVariable('col',np.int,('col'))
    day = np.arange(nday)
    sector = np.arange(nsector)
    row = np.arange(nrow)
    col = np.arange(ncol)
    fout.variables['day'][:] = day
    fout.variables['sector'][:] = sector
    fout.variables['row'][:] = row
    fout.variables['col'][:] = col
    fout.createVariable('MDA8SA', np.float32, ('day','sector','row','col'))
    fout.variables['MDA8SA'][:] = MDA8SA
    for num in range(len(Variables)):
        fout.createVariable(VariableNames[num], np.float32,('day','row','col'))
        fout.variables[VariableNames[num]][:] = Variables[num]
    fout.close()



##========================================================================##
##                            Input  Area                                 ##
##========================================================================##

  # 填写源解析的netCDF文件路径
O3SAfilePath = '/data3/liuhanqing/projects/postproc/2.combine/output/201811_OSAT/'    
  # 填写源解析的netCDF文件名                             
O3SAfileName = 'camx.YRD4km.Suzhou_fin.201811.OSAT.sa.grd01.ncf' 

  # 这一行填写未源解析netCDF文件路径
O3TOTfilePath = '/data3/liuhanqing/projects/postproc/2.combine/output/201811_OSAT/'   
  # 这一行填写未源解析netCDF文件名
O3TOTfilename = 'camx.YRD4km.Suzhou_fin.201811.OSAT.sa.grd01.ncf'   
  # 总臭氧小时浓度的变量名                      
O3TOTvarname = 'O3_R1R2R3'  

  # 输入你使用的模式(CMAQ/CAMx)
CTM = 'CAMx' 

  # 填写sector数、网格规格和源解析小时数
nsector = 14
nhour = 720                                                          
nrow = 216
ncol = 192


begin = datetime.datetime.now() 
print("START TIME: %s. It is estimated to take about 2-3 minutes, please be patient. " %begin)
if CTM == 'CAMx':
    # 将总臭氧浓度数据转化为八小时平均浓度。
    DA8TOT = treat_O3TOTdata(O3TOTfilePath,O3TOTfilename,O3TOTvarname,nhour,nrow,ncol)
    # 将臭氧源解析浓度数据转化为八小时平均浓度。
    DA8SA = treat_O3SAdata_CAMx(O3SAfilePath,O3SAfileName,nsector,nhour,nrow,ncol)
    # 根据总臭氧浓度的八小时平均值对臭氧源解析浓度数据进行解析，得到MDA8解析值。
    MDA8SA = cal_MDA8SA_CAMx(DA8TOT,DA8SA,nsector,nhour,nrow,ncol)

    # 标记需要获取的浓度贡献数据
    MDA8SA_NOx = np.sum(MDA8SA[:,0,:13,:,:],axis = 1)
    MDA8SA_VOC = np.sum(MDA8SA[:,1,:13,:,:],axis = 1)
    MDA8SA_IC = np.sum(MDA8SA[:,:,0,:,:],axis = 1)
    MDA8SA_BC = np.sum(MDA8SA[:,:,1,:,:],axis = 1)
    MDA8SA_PP = np.sum(MDA8SA[:,:,2,:,:],axis = 1)
    MDA8SA_IN = np.sum(MDA8SA[:,:,3,:,:],axis = 1)
    MDA8SA_PR = np.sum(MDA8SA[:,:,4,:,:],axis = 1)
    MDA8SA_DO = np.sum(MDA8SA[:,:,5,:,:],axis = 1)
    MDA8SA_TR = np.sum(MDA8SA[:,:,6,:,:],axis = 1)
    MDA8SA_DU = np.sum(MDA8SA[:,:,7,:,:],axis = 1)
    MDA8SA_OP = np.sum(MDA8SA[:,:,8,:,:],axis = 1)
    MDA8SA_AR = np.sum(MDA8SA[:,:,9,:,:],axis = 1)
    MDA8SA_MEGAN = np.sum(MDA8SA[:,:,10,:,:],axis = 1)
    MDA8SA_R2 = np.sum(MDA8SA[:,:,11,:,:],axis = 1)
    MDA8SA_R3 = np.sum(MDA8SA[:,:,12,:,:],axis = 1)
    MDA8SA_R2R3 = np.sum(MDA8SA[:,:,13,:,:],axis = 1)
    # 填写要标记的变量名(上下一致)
    VariableNames = ['MDA8SA_NOx','MDA8SA_VOC',
            'MDA8SA_IC','MDA8SA_BC',
            'MDA8SA_PP','MDA8SA_IN',
            'MDA8SA_PR','MDA8SA_DO',
            'MDA8SA_TR','MDA8SA_DU',
            'MDA8SA_OP','MDA8SA_AR',
            'MDA8SA_MEGAN','MDA8SA_R2',
            'MDA8SA_R3','MDA8SA_R2R3']
    Variables =  [MDA8SA_NOx,MDA8SA_VOC,
                    MDA8SA_IC,MDA8SA_BC,
                    MDA8SA_PP,MDA8SA_IN,
                    MDA8SA_PR,MDA8SA_DO,
                    MDA8SA_TR,MDA8SA_DU,
                    MDA8SA_OP,MDA8SA_AR,
                    MDA8SA_MEGAN,MDA8SA_R2,
                    MDA8SA_R3,MDA8SA_R2R3]
    # 将标记后的源解析结果写入netCDF文件中。
    save_MDA8SA_CAMx(MDA8SA,Variables,VariableNames,nhour,nrow,ncol)
    end = datetime.datetime.now()
    print ("END TIME: %s. The program was completed successfully! " %end)

elif CTM=='CMAQ':
        # 将总臭氧浓度数据转化为八小时平均浓度。
    DA8TOT = treat_O3TOTdata(O3TOTfilePath,O3TOTfilename,O3TOTvarname,nhour,nrow,ncol)
    # 将臭氧源解析浓度数据转化为八小时平均浓度。
    DA8SA = treat_O3SAdata_CMAQ(O3SAfilePath,O3SAfileName,nsector,nhour,nrow,ncol)
    # 根据总臭氧浓度的八小时平均值对臭氧源解析浓度数据进行解析，得到MDA8解析值。
    MDA8SA = cal_MDA8SA_CMAQ(DA8TOT,DA8SA,nsector,nhour,nrow,ncol)

    # 标记需要获取的浓度贡献数据
    MDA8SA_TOT = np.sum(MDA8SA[:,:,:,:],axis = 1)
    MDA8SA_BCO = np.sum(MDA8SA[:,0,:,:])
    MDA8SA_ICO = np.sum(MDA8SA[:,1,:,:])
    MDA8SA_SH_AR = np.sum(MDA8SA[:,2,:,:])
    MDA8SA_SH_BB = np.sum(MDA8SA[:,3,:,:])
    MDA8SA_SH_DU = np.sum(MDA8SA[:,4,:,:])
    MDA8SA_SH_EGU = np.sum(MDA8SA[:,5,:,:])
    MDA8SA_SH_IN = np.sum(MDA8SA[:,6,:,:])
    MDA8SA_SH_IN2 = np.sum(MDA8SA[:,7,:,:])
    MDA8SA_SH_MO = np.sum(MDA8SA[:,8,:,:])
    MDA8SA_SH_NA = np.sum(MDA8SA[:,9,:,:])
    MDA8SA_SH_RE = np.sum(MDA8SA[:,10,:,:])
    MDA8SA_SH_SO = np.sum(MDA8SA[:,11,:,:])
    MDA8SA_SH_ST = np.sum(MDA8SA[:,12,:,:])

    VariableNames = ['MDA8SA_TOT','MDA8SA_BCO',
            'MDA8SA_ICO','MDA8SA_SH_AR',
            'MDA8SA_SH_BB','MDA8SA_SH_DU',
            'MDA8SA_SH_EGU','MDA8SA_SH_IN',
            'MDA8SA_SH_IN2','MDA8SA_SH_MO',
            'MDA8SA_SH_NA','MDA8SA_SH_RE',
            'MDA8SA_SH_SO','MDA8SA_SH_ST']
    Variables =  [MDA8SA_TOT,MDA8SA_BCO,
            MDA8SA_ICO,MDA8SA_SH_AR,
            MDA8SA_SH_BB,MDA8SA_SH_DU,
            MDA8SA_SH_EGU,MDA8SA_SH_IN,
            MDA8SA_SH_IN2,MDA8SA_SH_MO,
            MDA8SA_SH_NA,MDA8SA_SH_RE,
            MDA8SA_SH_SO,MDA8SA_SH_ST]
    # 将标记后的源解析结果写入netCDF文件中。
    save_MDA8SA_CMAQ(MDA8SA,Variables,VariableNames,nhour,nrow,ncol)
    end = datetime.datetime.now()
    print ("END TIME: %s. The program was completed successfully! " %end)
else:
    print('You can only input CMAQ or CAMx to CTM, not %s. Program ends!'%CTM)




