import numpy as np
import netCDF4 as nc

####################################更新备注########################################
# 1.17     完成数据提取模块prep_var(path, filename)，
# 1.17     完成时间维度转化模块transVarTime(var_5d)，
# 1.17     完成读取MDA8时间索引模块get_MDA8TimeInd(path_ind,filename_ind,var),
# 1.17     完成源解析计算模块cal_MDA8sa(var1,var2) 
# 1.17     完成源解析结果写入模块mda8saAll_2nc(name,var1,var2,var3,varname1,varname2,varname3,list,var_list)
###################################################################################

def prep_var(path, filename):
    data = nc.Dataset(path + filename)
    varnames = list(data.variables.keys())  #read the nc variable name as a list
    varnames.remove('TFLAG')
    O3_5d = np.zeros((744,2,11,216,192),dtype = float)  
    for i in range(22):
        O3_5d[:,i%2,i//2,:,:] = data[varnames[i]] [23:,0,:,:]
    data.close()
    return O3_5d

def transVarTime(var_5d):
    var_6d = np.zeros((31,24,2,11,216,192), dtype=float)
    for day in range(31):
        for hour in range(24):
            var_6d[day,hour,:,:,:,:] = var_5d[24*day + hour,:,:,:,:]
    return var_6d

def get_MDA8TimeInd(path_ind,filename_ind,var):
    data = nc.Dataset(path_ind + filename_ind)
    mda8TimeInd = data[var][:]
    data.close()
    return mda8TimeInd

def cal_MDA8sa(var1,var2):
    mda8_sa_6d = np.zeros((31,24,2,11,216,192),dtype = float)
    for spec in range(2):
        for sector in range(11):
            mda8_sa_6d[:,:,spec,sector,:,:] = np.multiply(var1,var2[:,:,spec,sector,:,:])
    mda8_sa = np.zeros((31,2,11,216,192),dtype = float)
    for row in range(216):
        for col in range(192):
            for spec in range(2):
                for sector in range(11): 
                    for day in range(31):
                        for hour in range(24):
                            #print(hour)
                            #print(mda8_sa_6d[day,hour,spec,sector,row,col])
                            if mda8_sa_6d[day,hour,spec,sector,row,col] > 0.00:
                                # O3_sa_Ind = mda8_sa_6d.index(mda8_sa[day,hour,spec,sector,row,col])
                                #print(mda8_sa_6d[day,hour,spec,sector,row,col])
  
                                    
                                if (hour > 16)&(mda8_sa_6d[day,23,spec,sector,row,col] > 0):
                                    # print(mda8_sa_6d[day,:,spec,sector,row,col])
                                    # print('第%s行'%row)
                                    # print('第%s列'%col)
                                    # print('第%sday'%day)
                                    mda8_sa[day,spec,sector,row,col] = np.nanmean(np.hstack((mda8_sa_6d[day,hour:24,spec,sector,row,col],mda8_sa_6d[day + 1,0 : 8 - (24 - hour),spec,sector,row,col])))
                                    #print(hour)
                                else:
                                    mda8_sa[day,spec,sector,row,col] = np.nanmean(mda8_sa_6d[day,hour:hour + 8,spec,sector,row,col])
                                    #print(hour)
                                    
    return mda8_sa


def mda8sa_2nc(var,name,varname):
    fout = nc.Dataset(name,'w',format = 'NETCDF4')
    fout.createDimension('day',31)  
    # fout.createDimension('species',2)  
    # fout.createDimension('sector',11)  
    fout.createDimension('row',216)   
    fout.createDimension('col',192)
    fout.createVariable('day',np.int,('day'))  
    # fout.createVariable('species',np.int,('species'))
    # fout.createVariable('sector',np.int,('sector'))
    fout.createVariable('row',np.int,('row'))  ``
    fout.createVariable('col',np.int,('col'))
    day = np.arange(31)
    # species = np.arange(2)
    # sector = np.arange(11)
    row = np.arange(216)
    col = np.arange(192)
    fout.variables['day'][:] = day
    # fout.variables['species'][:] = species
    # fout.variables['sector'][:] = sector
    fout.variables['row'][:] = row
    fout.variables['col'][:] = col
    fout.createVariable(varname, np.float32, ('day','row','col'))
    # fout.createVariable(varname, np.float32, ('day','species','sector','row','col'))
    fout.variables[varname][:] = var
    fout.close()

def mda8saAll_2nc(name,var1,var2,var3,varname1,varname2,varname3,list,var_list):
    fout = nc.Dataset(name,'w',format = 'NETCDF4')
    fout.createDimension('day',31)  
    fout.createDimension('species',2)  
    fout.createDimension('sector',11)  
    fout.createDimension('row',216)   
    fout.createDimension('col',192)
    fout.createVariable('day',np.int,('day')) 
    fout.createVariable('species',np.int,('species'))
    fout.createVariable('sector',np.int,('sector'))  
    fout.createVariable('row',np.int,('row'))  
    fout.createVariable('col',np.int,('col'))
    day = np.arange(31)
    species = np.arange(2)
    sector = np.arange(11)
    row = np.arange(216)
    col = np.arange(192)
    fout.variables['day'][:] = day
    fout.variables['species'][:] = species
    fout.variables['sector'][:] = sector
    fout.variables['row'][:] = row
    fout.variables['col'][:] = col
    fout.createVariable(varname1, np.float32, ('day','species','sector','row','col'))
    fout.variables[varname1][:] = var1
    fout.createVariable(varname2, np.float32, ('day','row','col'))
    fout.variables[varname2][:] = var2
    fout.createVariable(varname3, np.float32, ('day','row','col'))
    fout.variables[varname3][:] = var3
    for j in range(11):
        fout.createVariable(list[j], np.float32, ('day','row','col'))
        fout.variables[list[j]][:] = var_list[j]
    fout.close()

###########################################
#              INPUT AREA                 #
###########################################

path = '/data3/liuhanqing/projects/postproc/2.combine/output/201807_OSAT/'
filename = 'camx.YRD4km.NOX_VOCs_new.201807.OSAT.sa.grd01.ncf'
path_ind = '/data/home/xuejin/code/python'
filename_ind = '/mda8_TimeInd.nc'

O3_5d = prep_var(path,filename)
O3_6d = transVarTime(O3_5d)
mda8_TimeInd = get_MDA8TimeInd(path_ind,filename_ind,'mda8TimeInd')
mda8_sa = cal_MDA8sa(mda8_TimeInd,O3_6d)

# mda8_sa = np.sum(O3_sa,axis = 1)/8
mda8_sa_nox = np.sum(mda8_sa[:,0,:,:,:],axis = 1)
mda8_sa_voc = np.sum(mda8_sa[:,1,:,:,:],axis = 1)
mda8_sa_IC = np.sum(mda8_sa[:,:,0,:,:],axis = 1)
mda8_sa_BC = np.sum(mda8_sa[:,:,1,:,:],axis = 1)
mda8_sa_PP = np.sum(mda8_sa[:,:,2,:,:],axis = 1)
mda8_sa_IN = np.sum(mda8_sa[:,:,3,:,:],axis = 1)
mda8_sa_PR = np.sum(mda8_sa[:,:,4,:,:],axis = 1)
mda8_sa_DO = np.sum(mda8_sa[:,:,5,:,:],axis = 1)
mda8_sa_TR = np.sum(mda8_sa[:,:,6,:,:],axis = 1)
mda8_sa_DU = np.sum(mda8_sa[:,:,7,:,:],axis = 1)
mda8_sa_OP = np.sum(mda8_sa[:,:,8,:,:],axis = 1)
mda8_sa_AR = np.sum(mda8_sa[:,:,9,:,:],axis = 1)
mda8_sa_MEGAN = np.sum(mda8_sa[:,:,10,:,:],axis = 1)
sector = ['mda8_sa_IC','mda8_sa_BC','mda8_sa_PP','mda8_sa_IN','mda8_sa_PR','mda8_sa_DO','mda8_sa_TR','mda8_sa_DU','mda8_sa_OP','mda8_sa_AR','mda8_sa_MEGAN']
mda8_sector =  [mda8_sa_IC,mda8_sa_BC,mda8_sa_PP,mda8_sa_IN,mda8_sa_PR,mda8_sa_DO,mda8_sa_TR,mda8_sa_DU,mda8_sa_OP,mda8_sa_AR,mda8_sa_MEGAN]
mda8saAll_2nc('mda8_sa_All.nc',mda8_sa,mda8_sa_nox,mda8_sa_voc,'mda8_sa','mda8_sa_nox','mda8_sa_voc',sector,mda8_sector)

mda8 = mda8_sa_nox + mda8_sa_voc
mda8sa_2nc(mda8,'mda8_new.nc','mda8')

#print(mda8.shape)









