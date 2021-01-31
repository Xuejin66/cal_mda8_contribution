import numpy as np
import netCDF4 as nc

####################################更新备注########################################
# 1.17     完成数据提取模块prep_var(path, filename)，
# 1.17     完成时间维度转化模块transVarTime(var_5d)，
# 1.17     完成读取MDA8时间索引模块get_MDA8TimeInd(path_ind,filename_ind,var),
# 1.17     完成源解析计算模块cal_MDA8sa(var1,var2) 
# 1.17     完成源解析结果写入模块mda8saAll_2nc(name,var1,var2,var3,varname1,varname2,varname3,list,var_list)
###################################################################################

###########################################
#           DEFINE FUNCTION               #
###########################################
def treat_unSAfile(unSAfilePath,unSAfilename,unSAvarname,nhour,nrow,ncol):
    unSAdata = nc.Dataset(unSAfilePath + unSAfilename)
    # varnames = list(data.variables.keys())  #read the nc variable name as a list
    unSApol = unSAdata[unSAvarname][:744,0,:,:]
    unSAdata.close()
    nday = nhour//24 - 1 
    unSApol_transTime = np.zeros((nday,32,nrow,ncol),dtype = float)
    for row in range(nrow):
        for col in range(ncol):
            for hour in range(nhour - 24):   #  只算30天的
                day = (hour + 1) // 24 
                unSApol_transTime[day,:,row,col] = unSApol[day*24:day*24 + 32,row,col]
    return unSApol_transTime
                


    # nday = nhour/24 -1
    # mda8 = np.zeros((nday,nrow,ncol),dtype = float)
    # mda8TimeInd = np.zeros((nday,nrow,ncol),dtype = int)
    # for col in range(192):
    #     for row in range(216):
    #         for hour in range(744):
    #             day = hour // 24 
    #             d8 = unSApol[day*24 : day*24 +8 , row, col]
    #             for hour in range(24):
    #                 d8 = var[hour:hour + 8, row, col]
    #                 da8 = np.nanmean(d8)
    #                 mda8_list.append(da8)    
    #             if len(mda8_list) == 24:                                                                                                                                                                                                                       
    #                 mda8[day,row,col] = max(mda8_list)
    #                 time_index = mda8_list.index(max(mda8_list))
    #                 mda8_TimeInd[day*24 + time_index:day*24 + time_index + 8,row,col] = 1
    #             print('时间%d'%time)
    #             print('row %d'%row)
    #             print('col %d'%col)
    # return mda8,mda8_TimeInd
    # return unSApol




###########################################
#              INPUT AREA                 #
###########################################

unSAfilePath = '/data3/liuhanqing/projects/postproc/2.combine/output/201807_OSAT/'    #这一行填写未源解析netCDF文件路径
unSAfilename = 'camx.YRD4km.R1R2R3.201807.OSAT.sa.grd01.ncf'                          #这一行填写未源解析netCDF文件名
unSAvarname = 'O3_R1R2R3'  
nhour = 744                                                           #这一行填写未源解析的总臭氧小时浓度
nrow = 216
ncol = 192








print(treat_unSAfile(unSAfilePath,unSAfilename,unSAvarname,nhour,nrow,ncol))
# # O3_4d = trans2var4d(O3_3d)
# mda8,mda8_TimeInd = cal_mda8(O3_3d)
# print(mda8_TimeInd[:,125,125])
# # O3_4d2nc(O3_4d,'O3_4d.nc')
# # mda8_2nc(mda8,'mda8.nc')
# mda8TimeInd2nc(mda8_TimeInd,'mda8_TimeInd.nc')

######