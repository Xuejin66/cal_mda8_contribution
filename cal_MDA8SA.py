import numpy as np
import netCDF4 as nc

##########################  FUNCTION EXPLAINATION  ################################
#
#      数据准备模块:         prep_SAdata(SAfilePath,SAfileName,nsector,nhour,nrow,ncol)，
#      MDA8时间索引模块:     get_MDA8TimeInd(path_ind,filename_ind,var),
#      源解析计算模块:       cal_MDA8sa(pol_6d,mda8Ind,nsector,nhour,nrow,ncol):
#      源解析结果写入模块:   mda8sa_2nc(mda8SA,nsector,nhour,nrow,ncol)
#
###################################################################################

def prep_SAdata(SAfilePath,SAfileName,nsector,nhour,nrow,ncol):
    data = nc.Dataset(SAfilePath + SAfileName)
    varnames = list(data.variables.keys())  #read the nc variable name as a list
    varnames.remove('TFLAG')
    pol_5d = np.zeros((nhour,2,nsector,nrow,ncol),dtype = float)  
    for i in range(2*nsector):
        pol_5d[:,i%2,i//2,:,:] = data[varnames[i]] [:nhour,0,:,:]
    data.close()
    nday = nhour // 24 - 1
    pol_6d = np.zeros((nday,31,2,nsector,nrow,ncol), dtype=float)
    for day in range(nday):
        # for hour in range(31):
        pol_6d[day,:,:,:,:,:] = pol_5d[24*day:24*day + 31,:,:,:,:]
    return pol_6d

def get_MDA8TimeInd(mda8IndPath,mda8IndFileName,mda8IndVarName):
    data = nc.Dataset(mda8IndPath + mda8IndFileName)
    mda8Ind = data[mda8IndVarName][:]
    data.close()
    return mda8Ind

def cal_MDA8sa(pol_6d,mda8Ind,nsector,nhour,nrow,ncol):
    nday = nhour // 24 - 1
    polSA_6d = np.zeros((nday,31,2,nsector,nrow,ncol),dtype = float)
    for spec in range(2):
        for sector in range(nsector):
            polSA_6d[:,:,spec,sector,:,:] = np.multiply(pol_6d[:,:,spec,sector,:,:],mda8Ind)
    mda8SA = np.sum(polSA_6d,axis = 1)/8
    return mda8SA

def mda8sa_2nc(mda8SA,nsector,nhour,nrow,ncol):
    nday = nhour // 24 - 1
    fout = nc.Dataset('mda8SA.nc','w',format = 'NETCDF4')
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
    fout.createVariable('mda8SA', np.float32, ('day','species','sector','row','col'))
    fout.variables['mda8SA'][:] = mda8SA
    fout.close()

def mda8_2nc(mda8,nhour,nrow,ncol):
    fout = nc.Dataset('mda8_SA.nc','w',format = 'NETCDF4')
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
    fout.variables['mda8'][:] = mda8
    fout.close()

def mda8SA_2nc(mda8SA,Variables,VariableNames,nhour,nrow,ncol):
    nday = nhour // 24 - 1
    fout = nc.Dataset('mdaSA.nc','w',format = 'NETCDF4')
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
    fout.createVariable('mda8SA', np.float32, ('day','species','sector','row','col'))
    fout.variables['mda8SA'][:] = mda8SA
    for num in range(len(Variables)):
        fout.createVariable(VariableNames[num], np.float32,('day','row','col'))
        fout.variables[VariableNames[num]][:] = Variables[num]
    fout.close()

###########################################
#              INPUT AREA                 #
###########################################
  # 填写源解析的netCDF文件路径
SAfilePath = '/data3/liuhanqing/projects/postproc/2.combine/output/201807_OSAT/'    
  # 填写源解析的netCDF文件名      
SAfileName = 'camx.YRD4km.NOX_VOCs_new.201807.OSAT.sa.grd01.ncf'                          
 
  # 填写上一个脚本计算出的mda8时间索引netCDF文件路径
mda8IndPath = '/data/home/xuejin/code/python/'                                                                 
  # 填写上一个脚本计算出的mda8时间索引netCDF文件名
mda8IndFileName = 'mda8Ind.nc'                                                            
  
  # 填写sector数、网格规格和源解析小时数
nsector = 11
nhour = 744                                                          
nrow = 216
ncol = 192


###########################################
#              MAIN PROGRAM               #
###########################################

  # 准备数据，将模型源解析得到的netCDF文件中的各类变量同化到一个六维列表中(list)。
pol_6d = prep_SAdata(SAfilePath,SAfileName,nsector,nhour,nrow,ncol)
  # 提取上个脚本中导出的时间索引数据。
mda8Ind = get_MDA8TimeInd(mda8IndPath,mda8IndFileName,'mda8Ind')
  # 计算解析结果，将时间索引数据×同化后的六维列表，按每日加和得到未标记的解析结果。
mda8SA = cal_MDA8sa(pol_6d,mda8Ind,nsector,nhour,nrow,ncol)
  # 将未标记的解析结果写入netCDF文件中。
mda8sa_2nc(mda8SA,nsector,nhour,nrow,ncol)
  # 分别按物种，人为源，天然源计算MDA8的贡献，得到标记后的源解析结果。
mda8SA_NOx = np.sum(mda8SA[:,0,:,:,:],axis = 1)
mda8SA_VOC = np.sum(mda8SA[:,1,:,:,:],axis = 1)
mda8SA_IC = np.sum(mda8SA[:,:,0,:,:],axis = 1)
mda8SA_BC = np.sum(mda8SA[:,:,1,:,:],axis = 1)
mda8SA_PP = np.sum(mda8SA[:,:,2,:,:],axis = 1)
mda8SA_IN = np.sum(mda8SA[:,:,3,:,:],axis = 1)
mda8SA_PR = np.sum(mda8SA[:,:,4,:,:],axis = 1)
mda8SA_DO = np.sum(mda8SA[:,:,5,:,:],axis = 1)
mda8SA_TR = np.sum(mda8SA[:,:,6,:,:],axis = 1)
mda8SA_DU = np.sum(mda8SA[:,:,7,:,:],axis = 1)
mda8SA_OP = np.sum(mda8SA[:,:,8,:,:],axis = 1)
mda8SA_AR = np.sum(mda8SA[:,:,9,:,:],axis = 1)
mda8SA_MEGAN = np.sum(mda8SA[:,:,10,:,:],axis = 1)
  # 填写要标记的变量名(上下一致，个数 = nsector + 2)
VariableNames = ['mda8SA_NOx','mda8SA_VOC',
          'mda8SA_IC','mda8SA_BC',
          'mda8SA_PP','mda8SA_IN',
          'mda8SA_PR','mda8SA_DO',
          'mda8SA_TR','mda8SA_DU',
          'mda8SA_OP','mda8SA_AR',
          'mda8SA_MEGAN']
Variables =  [mda8SA_NOx,mda8SA_VOC,
                mda8SA_IC,mda8SA_BC,
                mda8SA_PP,mda8SA_IN,
                mda8SA_PR,mda8SA_DO,
                mda8SA_TR,mda8SA_DU,
                mda8SA_OP,mda8SA_AR,
                mda8SA_MEGAN]
  # 将标记后的源解析结果写入netCDF文件中。
mda8SA_2nc(mda8SA,Variables,VariableNames,nhour,nrow,ncol)
mda8 = mda8SA_NOx + mda8SA_VOC
mda8_2nc(mda8,nhour,nrow,ncol)

