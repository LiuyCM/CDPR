#Created by Liu yifan 2023-05-26
#Connect me by email: liu.y.cm@m.titech.ac.jp
#Last modified time: 2023-06-30
#This file contains some functions that can be used to simulate

import numpy as np
from scipy.stats import norm
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation as R
import itertools
import math
import geometric_elements as ge

ASSEMBLING_ERROR_PARA = 0.5
AEP = ASSEMBLING_ERROR_PARA
DEFAULT_LINEAR_ERROR_RATIO = 0.02
DLER = DEFAULT_LINEAR_ERROR_RATIO
DEFAULT_CONFIDENCE_LEVEL = 0.90
DCL = DEFAULT_CONFIDENCE_LEVEL

ASSEMBLING_ERROR_BOUND = 50
AEB = ASSEMBLING_ERROR_BOUND
DEFAULT_LINEAR_ERROR_BOUND = 0.04
DLB = DEFAULT_LINEAR_ERROR_BOUND

3
#Function：产生正态分布的误差
#Input：interval：置信区间，confidence：置信度
def generate_norm_error(interval,confidence = None):
    if confidence is None:
        confidence = DCL
    z = norm.ppf(1-(1-confidence)/2)
    std_dev = interval/(2*z)
    error = np.random.normal(0, scale = std_dev)
    return error

#Function：产生线性误差
#Input：length：线性误差的长度，ratio：线性误差的比例，confidence：置信度
def generate_linear_error(length,ratio = None,confidence = None):
    if ratio is None:
        ratio = DLER
    if confidence is None:
        confidence = DCL
    interval = length*ratio
    error = generate_norm_error(interval,confidence)
    return error

#Function：产生装配误差
#Input：point_input：装配误差的生成点，confidence：置信度
def generate_point_assembling_error(point_input, interval = None,confidence = None,debug=None):
    if confidence is None:
        confidence = DCL
    if interval is None:
        interval = AEP
    point = point_input.copy()
    point.x += generate_linear_error(point.x,ratio = None,confidence = confidence) + generate_norm_error(interval,confidence = confidence)
    point.y += generate_linear_error(point.y,ratio = None,confidence = confidence) + generate_norm_error(interval,confidence = confidence)
    point.z += generate_linear_error(point.z,ratio = None,confidence = confidence) + generate_norm_error(interval,confidence = confidence)
    if debug == 'x':
        point.x += generate_norm_error(2*interval,confidence = confidence)
    return point

#Function: 模拟随机误差
#Input: point_input: 生成点，interval: 置信区间，confidence: 置信度
def generate_point_random_error(point_input,interval= None, confidence = None):
    if confidence is None:
        confidence = DCL
    if interval is None:
        interval = AEP
    point = point_input.copy()
    #point.debugcheck()
    point[0] += generate_norm_error(interval,confidence = confidence)
    point[1] += generate_norm_error(interval,confidence = confidence)
    point[2] += generate_norm_error(interval,confidence = confidence)
    #point.debugcheck()
    return point

#Function: 模拟编码器读取相对长度产生的误差
#Input: length: 相对长度，interval: 置信区间，confidence: 置信度
#因为我还没找到传感器的精度，暂且先认为传感器没有误差

#Below is the function to help calculate some parameters
#Function: 计算两点之间的距离
#Input: point1, point2: 两个点
def distance_between_two_points(point1,point2):
    distance = np.sqrt((point1[0]-point2[0])**2 + (point1[1]-point2[1])**2 + (point1[2]-point2[2])**2)
    return distance

#Function: 计算一个点所代表的向量的模
#Input: point: 一个点
def norm_of_point(point):
    norm = np.sqrt(point.x**2 + point.y**2 + point.z**2)
    return norm

#Function: 将一个点所代表的向量归一化
#Input: point: 一个点
def normalize_point(point):
    norm = norm_of_point(point)
    p = point.copy()
    p.x /= norm
    p.y /= norm
    p.z /= norm
    return p

#Function: 产生三角形内随机一个点的坐标，使用中心坐标法
#Input: point1, point2, point3: 三角形的三个顶点
def generate_point_in_triangle(point1,point2,point3):
    u = np.random.uniform(0,1)
    v = np.random.uniform(0,1)
    if u + v > 1:
        u = 1 - u
        v = 1 - v
    w = 1 - u - v
    point = point1.copy()
    point.x = point1.x*u + point2.x*v + point3.x*w
    point.y = point1.y*u + point2.y*v + point3.y*w
    point.z = point1.z*u + point2.z*v + point3.z*w
    return point

#Function: 根据传入的基点坐标，产生指定个必然在工作空间内的随机点
#Input: base_point: 基点坐标，num: 产生点的个数
#Output: a list of points
def generate_points_in_workspace(base_points,num):
    n = len(base_points)
    if n < 3 or num < 1:
        print("Error: the input is not correct")
        return
    triangle_list = []
    point_list = []
    for comb in itertools.combinations(base_points,3):
        triangle_list.append(comb)
    count = 0
    index = 0
    n = len(triangle_list)
    while count < num:
        triangle = triangle_list[index]
        point = generate_point_in_triangle(triangle[0],triangle[1],triangle[2])
        point_list.append(point)
        count += 1
        index = (index + 1) % n
    return point_list

def generate_lower_bound(x):
    x_l = []
    for i in x:
        i = (i - AEB - abs(i)*DLB)
        x_l.append(i)
    return x_l

def generate_upper_bound(x):
    x_u = []
    for i in x:
        i = (i + AEB + abs(i)*DLB)
        x_u.append(i)
    return x_u

def generate_bound(x):
    x_l = generate_lower_bound(x)
    x_u = generate_upper_bound(x)
    return (x_l,x_u)

#找出两个圆球的最小重力势能交点
def find_min_grav_2point_intersection(a, b, l1, l2):
    d = distance_between_two_points(a,b)

    if d > l1 + l2 or d <= abs(l1 - l2):
        # 两个圆球没有交点，返回到两个圆球表面距离之和最近的点
        return find_min_surface_distance_2point(a, b, l1, l2)
    else:  
        # 存在交点
        return find_triangle_point3d(a,b,l1,l2)

def find_triangle_point3d(a,b,l1,l2,minv = None):
    if minv is None:
        minv = -1
    d = distance_between_two_points(a,b)
    #theta is the angle between the line ab and the xy plane
    theta = math.atan2(b.z-a.z,math.sqrt((b.x-a.x)**2+(b.y-a.y)**2))
    #phi is the angle between the line ab and the yz plane
    phi = math.atan2(b.x-a.x,b.y-a.y)
    #lambda is the angle between the line ab and the line ac
    lamb = math.acos((l1**2+d**2-l2**2)/(2*l1*d))
    c = a.copy()
    c.x += l1*math.cos(theta+lamb*minv)*math.sin(phi)
    c.y += l1*math.cos(theta+lamb*minv)*math.cos(phi)
    c.z += l1*math.sin(theta+lamb*minv)
    return c

def find_min_surface_distance_2point(a, b, l1, l2):
    #ab代表从a指向b的单位矢量，ba同理
    ab = b - a
    ba = a - b
    ab = normalize_point(ab)
    ba = normalize_point(ba)
    # Return the midpoint
    c = a.copy()
    d = b.copy()
    c.x += ab.x * l1
    c.y += ab.y * l1
    c.z += ab.z * l1
    d.x += ba.x * l2
    d.y += ba.y * l2
    d.z += ba.z * l2
    c.x = (c.x + d.x) / 2
    c.y = (c.y + d.y) / 2
    c.z = (c.z + d.z) / 2
    return c

#计算正运动学的目标点，原理为把绳子当成弹簧，找到使得四根弹簧弹性势能最小的点，注意，下面两根
#绳子在小于初始长度时我们认为其弹性势能为零，以防出现反物理学的现象
def find_min_tension_intersection(base_points,cable_length,p = None):
    if(len(base_points.values()) !=  len(cable_length) and len(base_points.values()) != 4):
        print("Error: the FK input is not correct! The number of cable is not match the base points")
        return None
    if p == None:
        p = base_points[1].tolist()
    else:
        p = p.tolist()
    pointes = []
    for i in range(1,5):
        pointes.append(base_points[i].tolist())
    def objective(point):
        total_tension = 0
        for i in range(4):
            dis2center = distance_between_two_points(point,pointes[i])
            dis2surface = dis2center - cable_length[i]
            if dis2surface < 0:
                #为了防止出现下面两根线绷紧上面两根线松弛的反物理奇迹，我们假设上面两根线是弹簧
                if i == 1 or i == 2:
                    dis2surface = abs(dis2surface)
                else:
                    dis2surface = 0
            total_tension += dis2surface ** 2
        return total_tension
    res = minimize(objective, p, method='BFGS')
    ans = base_points[1].copy()
    ans.x = res.x[0]
    ans.y = res.x[1]
    ans.z = res.x[2]
    return ans

def workspace_projection(point,base_points):
    originP = base_points[1]
    def objective(x):
        newp = ge.Point3D(x[0],point.y+originP.y,point.z+originP.z)
        total_length = 0
        for p in base_points.values():
            total_length += distance_between_two_points(newp,p)
        return total_length
    res = minimize(objective, [point.x], method='BFGS')
    ans = ge.Point3D(res.x[0],point.y+originP.y,point.z+originP.z)
    return ans


def point_error(point,CDPR,mode,is2D):
    p = ge.Point3D(point[0],point[1],point[2])
    #pause
    cal_length = []
    if mode == "Before":
        length = CDPR.inverse_kinematic(IKmode = 1,MP_position=p,base_points= CDPR.cable_base_points_designed)
        for i in range(4):
            cal_length.append(length[i]-CDPR.cable_length_true[i]+CDPR.cable_length_init[i])
    elif mode == "After":
        p = workspace_projection(p,CDPR.cable_base_points_estimated)
        length = CDPR.inverse_kinematic(IKmode = 1,MP_position=p,base_points= CDPR.cable_base_points_estimated)
        for i in range(4):
            cal_length.append(length[i]-CDPR.cable_length_true[i]+CDPR.cable_length_estimated[i])
    elif mode == "After2D":
        p = workspace_projection(p,CDPR.cable_base_point2Ds_estimated)
        length = CDPR.inverse_kinematic(IKmode = 1,MP_position=p,base_points= CDPR.cable_base_point2Ds_estimated)
        for i in range(4):
            cal_length.append(length[i]-CDPR.cable_length_true[i]+CDPR.cable_length_estimated[i])
    elif mode == "True":
        p = workspace_projection(p,CDPR.cable_base_points)
        length = CDPR.inverse_kinematic(IKmode = 1,MP_position=p,base_points= CDPR.cable_base_points)
        for i in range(4):
            cal_length.append(length[i]-CDPR.cable_length_true[i]+CDPR.cable_length_true[i])
    else:
        print("Error: mode is not defined")
        return 99999
    resp = CDPR.forward_kinematic(FKmode=1,cable_length = cal_length)
    resp = resp - p
    if is2D != True:
        return resp.norm()
    else:
        return resp.norm(ignoreX = True)

def error_evaluation(CDPR,fig,ax,mode = None,is2D = None,quick=None):
    if mode == None:
        mode = "Before"
    if is2D == None:
        is2D = False
    start_y = 100
    end_y = 1300
    start_z = 100
    end_z = 3300
    if quick == False or quick == None:
        interval = 50
    else:
        interval = 100
    y = np.linspace(start_y,end_y,int((end_y - start_y+1)/interval))
    z = np.linspace(start_z,end_z,int((end_z - start_z+1)/interval))
    Y, Z = np.meshgrid(y,z)
    X = np.zeros_like(Y)
    Error = []
    totalnum = len(y)*len(z)
    counter = 0
    for xi, yi, zi in zip(np.ravel(X), np.ravel(Y), np.ravel(Z)):
        Error.append(point_error([xi, yi, zi],CDPR,mode,is2D = is2D))
        counter += 1
        if(counter % interval == 0):
            print(str(mode) + " Progress: " + "{:.3f}%".format(float(counter/totalnum)*100))
    Error = np.array(Error).reshape(X.shape)
    surface = ax.plot_surface(Y,Z,Error,cmap='viridis')
    #contour = ax.contourf(Y, Z, Error, levels=50, cmap='viridis')
    fig.colorbar(surface, label='Error')
    ax.set_xlabel('Y')
    ax.set_ylabel('Z')
    ax.set_zlabel('Error')
    pass

def get_rotation_matrix(euler,axis = None):
    #不传入轴则默认按照全局坐标旋转，传入euler角是一个三维数组分别表示xyz轴转角
    if axis is None:
        angle = [euler[0],euler[1],euler[2]]
        r = R.from_euler('xyz', angle, degrees=True)
    #传入轴则认为绕该轴旋转，且euler角是一个数，表示绕该轴转动的角度，角度制
    else:
        r_axis = np.array(axis)
        angle = np.deg2rad(euler)
        r = R.from_rotvec(angle * r_axis)
    return r.as_matrix()


def generate_length_sensor_error(length):
    #length是一个表示长度的标量，单位mm
    #返回一个标量，表示添加了长度传感器的误差后测得的长度值，单位mm
    #误差是一个高斯分布，均值为0，标准差为20微米每米，单位mm
    return length #+ length*np.random.normal(0,1e-5)

def generate_force_sensor_error(force):
    #force是一个表示力的标量，单位N
    #返回一个标量，表示添加了力传感器的误差后测得的力值，单位N
    #误差是一个高斯分布，均值为0，标准差为1%，单位N
    return force + force * np.random.normal(0,1e-2)