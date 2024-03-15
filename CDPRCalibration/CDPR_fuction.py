#Created by Liu yifan 2023-07-19
#Connect me by email: liu.y.cm@m.titech.ac.jp
#Last modified time: 2023-08-07
from distutils.log import debug
import numpy as np
import geometric_elements as ge
import time_utils as tu
from scipy.optimize import minimize
from scipy.optimize import least_squares

Average_force = 400


#逆运动学，输入CDPR和MP，根据输入条件返回末端执行器相对应的绳子连接点
def IK(CDPR,MP,fin = None,pulley = None):
    time_marker = tu.time_stamp()
    #同时有滑轮有鳍时的机构逆运动学
    if(fin == True and pulley == True):
        time_marker = tu.time_stamp()
        #time_marker1 = tu.time_stamp()
        P = MP.get_position()
        bi = MP.get_local_base_points()
        R = MP.get_rotation_matrix()
        RT = np.transpose(R)
        Ci = CDPR.get_cable_base_points()
        pBi = {}
        sigma = {}
        Ai = {}
        B_fin = {}
        #通过迭代鳍的转角求解逆运动学
        for i in range(1,5):
            def f(theta):
                theta = theta[0]
                #Step1: Assume a value of theta, and then calculate the position of the fin base point B_fin
                B_fin_local = np.array([0,MP.get_fin_length()*np.cos(theta),MP.get_fin_length()*np.sin(theta)])
                B_fin_local += bi[i].np()
                B_fin[i] = np.dot(R,B_fin_local) + P
                #Step2: use the position of B_fin to calculate the azimuth of i-th pulley sigma[i]
                #因为我们忽略了滑轮的方向误差默认其指向z轴正方，所以我们在这里不计算其旋转矩阵
                #Rp = CDPR.pulleys[i].get_rotation_matrix()
                #RpT = np.transpose(Rp)
                BCi = B_fin[i] - Ci[i].np()
                #pBi[i] = np.dot(RpT,BCi)
                pBi[i] = BCi
                sigma[i] = np.arctan2(pBi[i][1],pBi[i][0])
                CDPR.pulleys[i].set_tangent(np.degrees(sigma[i]))
                #Step3:Calculate the unit vector of i-th pulley
                u,w,k = CDPR.pulleys[i].get_unit_vector()
                #Step4:Create the equation about rotate angle of i-th pulley
                '''
                def g(x):
                    ni = np.cos(x)*u + np.sin(x)*k
                    lbat = B_fin[i] - Ci[i].np()-(CDPR.pulleys[i].radius*(u+ni))
                    return np.dot(lbat,ni)**2
                #Step5: use numerical method to calculate the rotate angle of i-th pulley
                res = minimize(g,0, method='BFGS',options={'disp': False})
                phi = res.x
                '''
                #solve a function about phi
                # a*np.cos(phi) + b*np.sin(phi) = c
                coe_c = CDPR.pulleys[i].radius
                coe_a = np.dot(BCi,u) - coe_c
                coe_b = np.dot(BCi,k)
                coe_R = np.sqrt(coe_a**2 + coe_b**2)
                if abs(coe_c) >= coe_R:
                    print("Warning! Mechanism interference!")
                    print("The distance between End-Effect fin and Pulley is less than radius")
                    raise("B_fin:"+ str(B_fin[i])+" C:" + str(Ci[i].np()))
                else:
                    #这里应用的是求一点对圆切点的方程，默认求的是右侧切点（也就是滑轮逆时针旋转放线），要根据i选择不同的解
                    alpha = np.arctan2(coe_b,coe_a)
                    if CDPR.pulleys[i].clock == True :
                        x1 = -np.arccos(coe_c / coe_R) + alpha
                    else:
                        x1 = np.arccos(coe_c / coe_R) + alpha
                    #最后别忘了舍去多余部分
                    phi = x1 % (2*np.pi)
                #Step6:Calculate the position of the connection point of i-th cable
                n = np.cos(phi)*u + np.sin(phi)*k
                Ai[i] = coe_c*(u + n) + Ci[i].np()
                #Step7: use the position of Ai to calculate the true position of B_fin
                A_local = Ai[i] - P
                A_local = np.dot(RT,A_local)
                true_theta = np.arctan2(A_local[2]-bi[i][2],A_local[1]-bi[i][1])
                #print(str(i)+" "+str(np.degrees(phi)))
                #In there we return the square of the error,because we do not use the least square method
                return (theta-true_theta)**2
            Ci_local = np.dot(RT,Ci[i].np()-P)
            theta = np.arctan2(Ci_local[2]-bi[i][2],Ci_local[1]-bi[i][1])
            #ini_theta = theta
            res = minimize(f,theta, method='BFGS',options={'disp': False})
            '''
            print("The "+str(i)+"-th fin has solved")
            print(res.message)
            print("The "+str(i)+"-th fin has rotated "+str(np.degrees(res.x))+" degrees")
            print("The "+str(i)+"-th fin has rotated from"+str(np.degrees(ini_theta))+" degrees")
            input("pause")
            '''
            #Debug information
            '''
            print("The "+str(i)+"-th fin has solved")
            print(res.message)
            print("The "+str(i)+"-th fin has rotated "+str(np.degrees(res.x))+" degrees")
            
            print("The "+str(i)+"-th cable base InverseKinematic takes " + tu.time_diff_str(time_marker))
            time_marker = tu.time_stamp()
            '''
        #print("InverseKinematic takes " + tu.time_diff_str(time_marker1))
        return B_fin, Ai
    #返回滑轮的转角，用于更新滑轮姿态
    elif(fin == None and pulley == True):
        #time_marker = tu.time_stamp()
        #time_marker1 = tu.time_stamp()
        P = MP.get_position()
        Bi = MP.get_fin_points()
        R = MP.get_rotation_matrix()
        Ci = CDPR.get_cable_base_points()
        pBi = {}
        sigma = {}
        angles = {}
        #print("InverseKinematic prepare takes " + tu.time_diff_str(time_marker))
        time_marker = tu.time_stamp()
        for i in range(1,5):
            #Step1:Calculate the azimuth of i-th pulley sigma[i]
            #Rp = CDPR.pulleys[i].get_rotation_matrix()
            #RpT = np.transpose(Rp)
            #pBi[i] = np.dot(RpT,Bi[i]-Ci[i].np())
            pBi[i] = Bi[i]-Ci[i].np()
            sigma[i] = np.degrees(np.arctan2(pBi[i][1],pBi[i][0]))
            CDPR.pulleys[i].set_rotation(sigma[i])
            #Step2:Calculate the unit vector of i-th pulley
            u,w,k = CDPR.pulleys[i].get_unit_vector()
            #Step3:Create the equation about rotate angle of i-th pulley
            '''
            def f(x):
                ni = np.cos(x)*u + np.sin(x)*k
                lbat = Bi[i] - Ci[i].np()-(CDPR.pulleys[i].radius*(u+ni))
                return np.dot(lbat,ni)**2
            #Step4:use numerical method to calculate the rotate angle of i-th pulley
            #phi is the rotate angle of i-th pulley
            #phi is our target
            phi = 0
            res = minimize(f,phi, method='BFGS',options={'disp': False})
            phi = res.x
            '''
            coe_c = CDPR.pulleys[i].radius
            coe_a = np.dot(pBi[i],u) - coe_c
            coe_b = np.dot(pBi[i],k)
            coe_R = np.sqrt(coe_a**2 + coe_b**2)
            if abs(coe_c) >= coe_R:
                print("Warning! Mechanism interference!")
                print("The distance between End-Effect fin and Pulley is less than radius")
                raise("B_fin:"+ str(B_fin[i])+" C:" + str(Ci[i].np()))
            else:
                #这里应用的是求一点对圆切点的方程，默认求的是右侧切点（也就是滑轮逆时针旋转放线），要根据i选择不同的解
                alpha = np.arctan2(coe_b,coe_a)
                if CDPR.pulleys[i].clock == True :
                    x1 = -np.arccos(coe_c / coe_R) + alpha
                else:
                    x1 = np.arccos(coe_c / coe_R) + alpha
                #最后别忘了舍去多余部分
                phi = x1 % (2*np.pi)
            angles[i]=phi
            '''
            #Step5:Calculate the position of the connection point of i-th cable
            n = np.cos(phi)*u + np.sin(phi)*k
            Ai[i] = CDPR.pulleys[i].radius*(u + n) + Ci[i].np()
            Ai[i] = ge.Point3D(Ai[i][0],Ai[i][1],Ai[i][2])
            
            #Debug information
            print("The "+str(i)+"-th pulley has solved")
            print(res.message)
            print("The "+str(i)+"-th pulley has rotated "+str(np.degrees(phi))+" degrees")
            print("The "+str(i)+"-th pulley has " + str(res.fun))
            print("The "+str(i)+"-th cable base point on EE is "+str(Bi[i]))
            '''
            #print("The "+str(i)+"-th cable base InverseKinematic takes " + tu.time_diff_str(time_marker))
            #time_marker = tu.time_stamp()
        #print("InverseKinematic ignore fin takes " + tu.time_diff_str(time_marker1))
        return angles
    #如果只有末端执行器，没有滑轮
    #返回鳍和绳子的链接点Bfin
    elif(fin == True and pulley == None):
        A = CDPR.get_cable_base_points()
        P = MP.get_position()
        R = MP.get_rotation_matrix()
        RT = np.transpose(R)
        A_local = {}
        B_local = MP.get_local_base_points()
        B_fin_local = {}
        B_fin = {}
        for i in range(1,5):
            A_local[i] = np.array([A[i].x - P[0],A[i].y - P[1],A[i].z - P[2]])
            A_local[i] = np.dot(RT,A_local[i])
            theta = np.arctan2(A_local[i][1]-B_local[i].y,A_local[i][2]-B_local[i].z)
            B_fin_local[i] = ge.Point3D(0,MP.get_fin_length()*np.sin(theta),MP.get_fin_length()*np.cos(theta))
            B_fin_local[i] += B_local[i]
            B_fin[i] = np.dot(R,B_fin_local[i].np())
            B_fin[i] += P
            B_fin[i] = ge.Point3D(B_fin[i][0],B_fin[i][1],B_fin[i][2])
        print("InverseKinematic ignore pulley takes " + tu.time_diff_str(time_marker))
        return B_fin
        
    else:
        print("WIP,please wait for update")
        raise NotImplementedError

#Inverse Statics
def IS(CDPR,MP):
    #time_marker = tu.time_stamp()
    #W is the force and moment matrix of the end effector in the world frame caused by gravity  
    #W = np.array([0,0,-MP.get_mass()*9.8,0,0,0,Average_force])
    W = np.array([0,MP.get_mass()*9.80665,0,Average_force])
    W = np.transpose(W)
    ki = {}
    bfi = MP.get_fin_points()
    Ai = CDPR.get_cable_attach_points()
    P = MP.get_position()
    #u is the unit vector of the cable
    ui = {}
    for i in range(1,5):
        ui[i] = Ai[i]-bfi[i]
        ui[i] = ui[i]/np.linalg.norm(ui[i])
        ki[i] = np.cross(bfi[i]-P,ui[i])
        #print("debugflag:",np.linalg.norm(ki[i]))
    #3x4 coefficient matrix 用于计算残差的矩阵
    B = np.array([[ui[1][0],ui[2][0],ui[3][0],ui[4][0]],
                  [ki[1][1],ki[2][1],ki[3][1],ki[4][1]],
                  [ki[1][2],ki[2][2],ki[3][2],ki[4][2]]])
    #4x4 coefficient matrix 用于计算力矩的矩阵
    A = np.array([[ui[1][1],ui[2][1],ui[3][1],ui[4][1]],
                  [ui[1][2],ui[2][2],ui[3][2],ui[4][2]],
                  [ki[1][0],ki[2][0],ki[3][0],ki[4][0]],
                  [1,1,1,1]])
    try:
        # 尝试计算逆矩阵
        inverse_matrix = np.linalg.inv(A)
    except np.linalg.LinAlgError:
        # 如果出现奇异点或不可逆情况，计算伪逆矩阵
        inverse_matrix = np.linalg.pinv(A)
        print("Pseudo-inverse Matrix(ignore moment)")
    force = np.dot(inverse_matrix,W)
    force = np.transpose(force)
    for f in force:
        if(f<0):
            print(force)
            print("Warning:Force is negative")
            print(A)
            print(inverse_matrix)
            print(W)
            input()
    residual = np.dot(B,force)
    #print(B)
    #print("InverseStatic takes " + tu.time_diff_str(time_marker))
    #返回力和残差
    return (force,residual)

#Workspace
def WS(CDPR,MP):
    time_markerWS = tu.time_stamp()
    pose = MP.get_pose()
    y = pose[1]
    z = pose[2]
    phi = pose[3]
    #输入MP的位置，返回此时的残差
    def EEposition(x):
        MP.set_pose([x[0],y,z,phi,x[1],x[2]])
        #print([x[0],y,z,phi,x[1],x[2]])
        #MP.update_base_points()
        B_fin,Ai = IK(CDPR,MP,pulley=True,fin=True)
        MP.set_fin_points(B_fin)
        CDPR.set_cable_attach_points(Ai)
        res = IS(CDPR,MP)
        residual = res[1]
        #r0单位是N， r1和r2单位是Nmm，所以需要进行1000倍的缩放
        residual[0] = residual[0]
        residual[1] = residual[1]/1e3
        residual[2] = residual[2]/1e3
        #print(residual)
        return residual
    # 初始优化变量 分别为x轴位置，绕y轴转角和绕z轴转角（角度制）
    x0 = np.array([0,0,0])
    # 定义优化变量的范围
    bnds = ([-200,-45,-45],[200,45,45])
    #result = minimize(EEposition,x0, method='BFGS',bounds=bnds,options={'disp': True})
    result = least_squares(EEposition, x0,ftol=1e-6, xtol=1e-6, gtol=1e-6,bounds=bnds,verbose=False)
    #print(result)
    print("Find a projection in workspace takes " + tu.time_diff_str(time_markerWS))
    return result

#Inverse calculation, used for calculate the residual
def IC(CDPR,MP,pose):
    #time_marker = tu.time_stamp()
    MP.set_pose(pose)
    P = pose[:3]
    bi = MP.get_local_base_points()
    R = MP.get_rotation_matrix()
    RT = np.transpose(R)
    Ci = CDPR.get_cable_base_points()
    pBi = {}
    sigma = {}
    Ai = {}
    B_fin = {}
    angles = {}
    tphi = [0]
    #通过迭代鳍的转角求解逆运动学
    for i in range(1,5):
        def f(theta):
            theta = theta[0]
            B_fin_local = np.array([0,MP.get_fin_length()*np.cos(theta),MP.get_fin_length()*np.sin(theta)])
            B_fin_local += bi[i].np()
            B_fin[i] = np.dot(R,B_fin_local) + P
            #Rp = CDPR.pulleys[i].get_rotation_matrix()
            #RpT = np.transpose(Rp)
            BCi = B_fin[i] - Ci[i].np()
            #pBi[i] = np.dot(RpT,BCi)
            pBi[i] = BCi
            sigma[i] = np.arctan2(pBi[i][1],pBi[i][0])
            CDPR.pulleys[i].set_tangent(np.degrees(sigma[i]))
            u,w,k = CDPR.pulleys[i].get_unit_vector()
            coe_c = CDPR.pulleys[i].radius
            coe_a = np.dot(BCi,u) - coe_c
            coe_b = np.dot(BCi,k)
            coe_R = np.sqrt(coe_a**2 + coe_b**2)
            if abs(coe_c) > coe_R:
                print("Warning! Mechanism interference!")
                print("The distance between End-Effect fin and Pulley is less than radius")
                raise("B_fin:"+ str(B_fin[i])+" C:" + str(Ci[i].np()))
            else:
                alpha = np.arctan2(coe_b,coe_a)
                if CDPR.pulleys[i].clock == True :
                    x1 = -np.arccos(coe_c / coe_R) + alpha
                else:
                    x1 = np.arccos(coe_c / coe_R) + alpha
                phi = x1 % (2*np.pi)
            #Step6:Calculate the position of the connection point of i-th cable
            n = np.cos(phi)*u + np.sin(phi)*k
            Ai[i] = coe_c*(u + n) + Ci[i].np()
            #Step7: use the position of Ai to calculate the true position of B_fin
            A_local = Ai[i] - P
            A_local = np.dot(RT,A_local)
            true_theta = np.arctan2(A_local[2]-bi[i][2],A_local[1]-bi[i][1])
            #print(str(i)+" "+str(np.degrees(phi)))
            tphi[0] = phi
            return (theta-true_theta)**2
        Ci_local = np.dot(RT,Ci[i].np()-P)
        theta = np.arctan2(Ci_local[2]-bi[i][2],Ci_local[1]-bi[i][1])
        res = minimize(f,theta,method='BFGS',options={'disp': False})
        angles[i] = tphi[0]
    MP.set_fin_points(B_fin)
    CDPR.set_cable_attach_points(Ai)
    ui = {}
    ki = {}
    for i in range(1,5):
        CDPR.pulleys[i].update(attachpoint = (Ai[i]))
        #print(np.degrees(abs(np.pi-abs(angles[i]))))
        CDPR.pulleys[i].attach_angle = abs(np.pi-abs(angles[i]))
        ui[i] = Ai[i] - B_fin[i]
        ui[i] = ui[i]/np.linalg.norm(ui[i])
        ki[i] = np.cross(B_fin[i]-P,ui[i])

    W = np.array([0,MP.get_mass()*9.80665,0,Average_force])
    W = np.transpose(W)
    #3x4 coefficient matrix 用于计算残差的矩阵
    B = np.array([[ui[1][0],ui[2][0],ui[3][0],ui[4][0]],
                  [ki[1][1],ki[2][1],ki[3][1],ki[4][1]],
                  [ki[1][2],ki[2][2],ki[3][2],ki[4][2]]])
    #4x4 coefficient matrix 用于计算力矩的矩阵
    A = np.array([[ui[1][1],ui[2][1],ui[3][1],ui[4][1]],
                  [ui[1][2],ui[2][2],ui[3][2],ui[4][2]],
                  [ki[1][0],ki[2][0],ki[3][0],ki[4][0]],
                  [1,1,1,1]])
    try:
        # 尝试计算逆矩阵
        inverse_matrix = np.linalg.inv(A)
    except np.linalg.LinAlgError:
        # 如果出现奇异点或不可逆情况，计算伪逆矩阵
        inverse_matrix = np.linalg.pinv(A)
        print("Pseudo-inverse Matrix(ignore moment)")

    forces = np.dot(inverse_matrix,W)
    forces = np.transpose(forces)
    residual = np.dot(B,forces)

    CDPR.update_force_sensors(forces=forces) 
    encoders = CDPR.read_encoders()
    #print("Inverse Calculation takes " + tu.time_diff_str(time_marker))
    return encoders,forces,residual

