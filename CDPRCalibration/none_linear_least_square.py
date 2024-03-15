#Created by Liu yifan 2023-05-27
#Connect me by email: liu.y.cm@titech.ac.jp
#Last modified time: 2023-08-05
#This file contains some functions that can be used to solve none linear least square problem about CDPR
#该文件为旧版本实验用，现在已经无法运行，保留供参考用
#Warning: All the functions in this file can only be used in the simplified CDPR model

#import numpy as np
import draw_utils as du
import geometric_elements as ge
from scipy.optimize import least_squares
import simulation_utils as su
import post_process as pp
import time_utils as tu



def from_list_to_vaiables(x,issimplified = False):
    cable_base_points = {}
    for i in range(4):
        cable_base_points[i] = ge.Point3D(x[i*3],x[i*3+1],x[i*3+2])
    cable_init_length = {}
    for i in range(4):
        cable_init_length[i] = x[12+i]
    if issimplified == False:
        sampling_points = {}
        for i in range((len(x)-16)//3):
            sampling_points[i] = ge.Point3D(x[16+i*3],x[16+i*3+1],x[16+i*3+2])
        return cable_base_points,cable_init_length,sampling_points
    else:
        return cable_base_points,cable_init_length

#sampling function
#input:the CDPR mechanism which need to be sampled and the sampling times(optional)
#output:a list contains the sampling result(the position of the end effector)
def generate_sampling(CDPR,num = None):
    if num == None:
        num = 50
    sampling = CDPR.random_sampling(num)
    return sampling

def generate_sampling_with_error(sampling_record):
    sampling_with_error = {}
    for i in range(len(sampling_record)):
        sampling_with_error[i] = su.generate_point_random_error(sampling_record[i])
    return sampling_with_error


#target function
#input:the position of  the end effector and the index of the cable
#output:the length of the cable
def nonlinear_function(CDPR,index,base_points = None,MP_position = None):
    cable_length = CDPR.inverse_kinematic(IKmode = 2,cable_index = index,base_points = base_points,MP_position = MP_position)
    return cable_length


#residual function
def residual(params,CDPR,y):
    cable_base_points,cable_init_length,sampling_points = from_list_to_vaiables(params,issimplified = False)
    resi = []
    for i in sampling_points.keys():
        for j in range(4):
            length_ij = CDPR.inverse_kinematic(IKmode = 2,cable_index = j,base_points = cable_base_points,MP_position = sampling_points[i])
            relative_length_ij = length_ij - cable_init_length[j]
            resi_ij =  relative_length_ij - y[i][j+1]
            resi.append(resi_ij)
    return resi


def residual_2(params,CDPR,y,sampling_points):
    cable_base_points,cable_init_length = from_list_to_vaiables(params,issimplified = True)
    resi = []
    for i in sampling_points.keys():
        for j in range(4):
            length_ij = CDPR.inverse_kinematic(IKmode = 2,cable_index = j,base_points = cable_base_points,MP_position = sampling_points[i])
            relative_length_ij = length_ij - cable_init_length[j]
            resi_ij =  relative_length_ij - y[i][j+1]
            resi.append(resi_ij)
    return resi

#create the x vector(unknown values)
def createx(CDPR,sam_num,sampling_record = None,length_init_true = None):
    x = []
    
    if sampling_record == None:
        for point in CDPR.cable_base_points_designed.values():
            x.append(point.x)
            x.append(point.y)
            x.append(point.z)
        for index in CDPR.cable_base_points_designed.keys():
            x.append(nonlinear_function(CDPR,index,base_points = CDPR.cable_base_points_designed,MP_position = CDPR.move_platfrom_position_initial))
        for i in range(sam_num):
            x.append(CDPR.move_platfrom_position_initial.x)
            x.append(CDPR.move_platfrom_position_initial.y)
            x.append(CDPR.move_platfrom_position_initial.z)
    elif length_init_true != None:
        for point in CDPR.cable_base_points.values():
            x.append(point.x)
            x.append(point.y)
            x.append(point.z)
        for index in CDPR.cable_base_points.keys():
            x.append(length_init_true[index])
        for i in range(sam_num):
            x.append(sampling_record[i].x)
            x.append(sampling_record[i].y)
            x.append(sampling_record[i].z)
    elif length_init_true == None:
        for point in CDPR.cable_base_points_designed.values():
            x.append(point.x)
            x.append(point.y)
            x.append(point.z)
        for index in CDPR.cable_base_points.keys():
            x.append(nonlinear_function(CDPR,index,base_points = CDPR.cable_base_points_designed,MP_position = CDPR.move_platfrom_position_initial))
        for i in range(sam_num):
            x.append(sampling_record[i].x)
            x.append(sampling_record[i].y)
            x.append(sampling_record[i].z)
    return x

#calibration function
#input:the CDPR mechanism which need to be calibrated and the sampling times(optional)
def NLLS_calibration(CDPR,ax,sam_num = None,do_post_process = None):
    
    time_marker = tu.time_stamp()
    length_init_true = {}
    for index in CDPR.cable_base_points.keys():
        length_init_true[index] = nonlinear_function(CDPR,index)
    #print(length_init_true)

    #生成采样点并将其储存在sampling_list内
    sampling_list = generate_sampling(CDPR,sam_num)
    #只绘制最多32个采样点
    if sam_num > 32:
        du.draw_points(sampling_list[:32], ax,color='k')
    else:
        du.draw_points(sampling_list, ax,color='k')

    sampling_record = {}
    for i in range(sam_num):
        sampling_record[i] = sampling_list[i]
    #生成带有误差的采样点坐标将其作为优化的起始值
    sampling_record_with_error = generate_sampling_with_error(sampling_record)
    for index, point in sampling_record_with_error.items():
        sampling_record_with_error[index] = ge.Point3D(0,point.y,point.z)


    x_true =  createx(CDPR,sam_num,sampling_record,length_init_true)
    
    
    # x means the variable to be optimized
    x = createx(CDPR,sam_num,sampling_record_with_error)
    # y means the variable we have known
    y={}
    for index in sampling_record.keys():
        y[index]=CDPR.simu_encoder(length_init_true,sampling_record[index])

    #resi = residual_2(x[:16],CDPR,y,sampling_record)
    #print(resi)
    #resi = residual_2(x_true[:16],CDPR,y,sampling_record)
    #print(resi)
    #result = least_squares(residual_2, x[:16], args=(CDPR, y,sampling_record))
    result = least_squares(residual, x, args=(CDPR, y),bounds = su.generate_bound(x))
    timeused = tu.time_diff(time_marker)
    postProcess = pp.post_processing(result=result,x = x,y = y,x_true = x_true,timeused = timeused/1000)
    if do_post_process == None:
        return result.x
    else:
        return postProcess
    
    
