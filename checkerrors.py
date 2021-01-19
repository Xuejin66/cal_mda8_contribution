import numpy as np
import netCDF4 as nc

def checkerrors_3d(path,filename,var):
    data = nc.Dataset(path + filename)
    checkvar = data[var][:,125,125]
    data.close()
    print(checkvar)

def checkerrors_4d(path,filename,var):
    data = nc.Dataset(path + filename)
    checkvar = data[var][23:,0,:,:]
    data.close()
    return checkvar

def checkerrors_4D(path,filename,var):
    data = nc.Dataset(path + filename)
    checkvar = data[var][:,:,:,:]
    data.close()
    return checkvar

path = "/data/home/xuejin/code/python/"
print('mda8，数据来源于camx.YRD4km.R1R2R3.201807.OSAT.sa.grd01.ncf')
checkerrors_3d(path,'mda8.nc','mda8')
print('源解析加和后的mda8,即下面两个的加和，数据为camx.YRD4km.NOX_VOCs_new.201807.OSAT.sa.grd01.ncf经过一系列数据处理')
checkerrors_3d(path,'mda8_new.nc','mda8')
print('mda8按nox源解析，数据为camx.YRD4km.NOX_VOCs_new.201807.OSAT.sa.grd01.ncf经过一系列数据处理')
checkerrors_3d(path,'mda8_sa_nox.nc','mda8_sa_nox')
print('mda8按voc源解析，数据为camx.YRD4km.NOX_VOCs_new.201807.OSAT.sa.grd01.ncf经过一系列数据处理')
checkerrors_3d(path,'mda8_sa_voc.nc','mda8_sa_voc')
varname = ['O3_IC_NOX','O3_BC_NOX','O3_R1G1_NOX','O3_R1G2_NOX','O3_R1G3_NOX','O3_R1G4_NOX','O3_R1G5_NOX','O3_R1G6_NOX','O3_R1G7_NOX','O3_R1G8_NOX','O3_R1G9_NOX'] 
O3_nox = np.zeros((744,216,192),dtype = float)
for spec in varname:
    temp = checkerrors_4d(path,'camx.YRD4km.NOX_VOCs_new.201807.OSAT.sa.grd01.ncf',spec)
    O3_nox = O3_nox + temp
print('7月7日 O3按nox源解析，数据为camx.YRD4km.NOX_VOCs_new.201807.OSAT.sa.grd01.ncf各变量加和')
print(O3_nox[144:168,125,125])

mda8TimeInd = checkerrors_4D(path,'mda8_TimeInd.nc','mda8TimeInd')
print('7月7日 O3源解析时间索引，数据来源于camx.YRD4km.R1R2R3.201807.OSAT.sa.grd01.ncf')
print(mda8TimeInd[6,:,125,125])





#O3_IC_NOX = checkerrors_4d(path,'camx.YRD4km.NOX_VOCs_new.201807.OSAT.sa.grd01.ncf','O3_IC_NOX')
#O3_BC_NOX = checkerrors_4d(path,'camx.YRD4km.NOX_VOCs_new.201807.OSAT.sa.grd01.ncf','O3_BC_NOX')
#O3_R1G1_NOX = checkerrors_4d(path,'camx.YRD4km.NOX_VOCs_new.201807.OSAT.sa.grd01.ncf','O3_R1G1_NOX')
#O3_R1G2_NOX = checkerrors_4d(path,'camx.YRD4km.NOX_VOCs_new.201807.OSAT.sa.grd01.ncf','O3_R1G2_NOX')
#O3_R1G3_NOX = checkerrors_4d(path,'camx.YRD4km.NOX_VOCs_new.201807.OSAT.sa.grd01.ncf','O3_R1G3_NOX')
#O3_R1G4_NOX = checkerrors_4d(path,'camx.YRD4km.NOX_VOCs_new.201807.OSAT.sa.grd01.ncf','O3_R1G4_NOX')
#O3_R1G5_NOX = checkerrors_4d(path,'camx.YRD4km.NOX_VOCs_new.201807.OSAT.sa.grd01.ncf','O3_R1G5_NOX')
#O3_R1G6_NOX = checkerrors_4d(path,'camx.YRD4km.NOX_VOCs_new.201807.OSAT.sa.grd01.ncf','O3_R1G6_NOX')
#O3_R1G7_NOX = checkerrors_4d(path,'camx.YRD4km.NOX_VOCs_new.201807.OSAT.sa.grd01.ncf','O3_R1G7_NOX')
#O3_R1G8_NOX = checkerrors_4d(path,'camx.YRD4km.NOX_VOCs_new.201807.OSAT.sa.grd01.ncf','O3_R1G8_NOX')
#O3_R1G9_NOX = checkerrors_4d(path,'camx.YRD4km.NOX_VOCs_new.201807.OSAT.sa.grd01.ncf','O3_R1G9_NOX')
#O3_nox = ['O3_IC_NOX','O3_BC_NOX','O3_R1G1_NOX','O3_R1G2_NOX','O3_R1G3_NOX','O3_R1G4_NOX','O3_R1G5_NOX','O3_R1G6_NOX','O3_R1G7_NOX','O3_R1G8_NOX','O3_R1G9_NOX'] 
#O3_nox = 





