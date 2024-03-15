#Created by Liu yifan 2023-11-12
#Connect me by email: liu.y.cm@m.titech.ac.jp
#Last modified time: 2023-11-13
import CDPR_fuction as CF
import numpy as np
import torch
import torchgeometry as tgm

#Inverse calculation, used for calculate the jacobian matrix (length part)
def inv_jac_length(CDPR,MP,pose):
    MP.set_pose(pose)
    P = pose[:3]
    B_fin,Ai = CF.IK(CDPR,MP,fin = True,pulley = True)
    b_fin = {}
    u = {}
    #this first ont is a 4*12 matrix, each row is the jacobian matrix of cable base point
    #the second one is a 4*4 matrix, each row is the jacobian matrix of initial cable length, this is a minus identity matrix
    #the third one is a 4*6 matrix, each row is the jacobian matrix of pose
    jac_1 = np.zeros((4,12))
    jac_2 = -1 * np.identity(4)
    jac_3 = np.zeros((4,6))
    for i in range(4):
        #bfin is the vector from the base point of MP to the end point of the fin
        #u is the vector from the end point of the fin to the connection point of the cable and pulley
        b_fin[i] = B_fin[i+1] - P
        u[i] = Ai[i+1] - B_fin[i+1]
        print("debug")
        print(Ai[i+1])
        print(B_fin[i+1])
        print(u[i])
        u[i] = u[i] / np.linalg.norm(u[i])
        print(u[i])
        jac_1[i][i*3] = u[i][0]
        jac_1[i][i*3+1] = u[i][1]
        jac_1[i][i*3+2] = u[i][2]
        radius = np.linalg.norm(b_fin[i])
        b_fin_unit = b_fin[i] / radius
        #this coefficient is used to convert the length to angle in degree
        coe = np.pi * radius/180
        omega = np.concatenate((u[i],coe * np.cross(b_fin_unit,u[i])))
        jac_3[i,:] = -1*omega
    return jac_1,jac_2,jac_3


def degrees_to_radians(degrees):
    return degrees * np.pi / 180

def radians_to_degrees(radians):
    return radians * 180 / np.pi

def inv_jac_force_2D(pose):
    #W is the wrench position, which is a 3*1 vector
    W = []
    W.append(torch.tensor([0,0,0]))
    W.append(torch.tensor([0,0,3640]))
    W.append(torch.tensor([0,1630,3640]))
    W.append(torch.tensor([0,1630,0]))
    EEwidth = 86
    EEheight = 49
    fin_length = 287
    P = torch.tensor(pose[:3])
    orientation = torch.tensor(pose[3:])
    orientation = degrees_to_radians(orientation)
    B = []
    B.append(np.array([P[0],P[1]-EEwidth/2,P[2]+EEheight/2]))
    B.append(np.array([P[0],P[1]-EEwidth/2,P[2]-EEheight/2]))
    B.append(np.array([P[0],P[1]+EEwidth/2,P[2]-EEheight/2]))
    B.append(np.array([P[0],P[1]+EEwidth/2,P[2]+EEheight/2]))
    RotationMatrix = tgm.euler_angles_to_matrix(orientation, "XYZ")

    R = MP.get_rotation_matrix()
    Ci = CDPR.get_cable_base_points()
    pBi = {}
    sigma = {}
    angles = {}
    #print("InverseKinematic prepare takes " + tu.time_diff_str(time_marker))
    time_marker = tu.time_stamp()
    for i in range(1,5):
        #Step1:Calculate the azimuth of i-th pulley sigma[i]
        Rp = CDPR.pulleys[i].get_rotation_matrix()
        RpT = np.transpose(Rp)
        pBi[i] = np.dot(RpT,Bi[i]-Ci[i].np())
        sigma[i] = np.degrees(np.arctan2(pBi[i][1],pBi[i][0]))
        CDPR.pulleys[i].set_rotation(sigma[i])
        #Step2:Calculate the unit vector of i-th pulley
        u,w,k = CDPR.pulleys[i].get_unit_vector()
        #Step3:Create the equation about rotate angle of i-th pulley
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
        angles[i]=phi
            
        #Step5:Calculate the position of the connection point of i-th cable
        n = np.cos(phi)*u + np.sin(phi)*k
        Ai[i] = CDPR.pulleys[i].radius*(u + n) + Ci[i].np()
        Ai[i] = ge.Point3D(Ai[i][0],Ai[i][1],Ai[i][2])
    B_fin,Ai = CF.IK(CDPR,MP,fin = True,pulley = True)
    P = pose[:3]
    W = np.transpose(W)
    ki = {}
    #u is the unit vector of the cable
    ui = {}
    for i in range(1,5):
        ui[i] = Ai[i]-B_fin[i]
        ui[i] = ui[i]/np.linalg.norm(ui[i])
        ki[i] = np.cross(B_fin[i]-P,ui[i])
    #4x4 coefficient matrix 
    A = np.array([[ui[1][1],ui[2][1],ui[3][1],ui[4][1]],
                  [ui[1][2],ui[2][2],ui[3][2],ui[4][2]],
                  [ki[1][0],ki[2][0],ki[3][0],ki[4][0]],
                  [1,1,1,1]])
    inverse_matrix = np.linalg.inv(A)
    force = np.dot(inverse_matrix,W)
    force = np.transpose(force)
    return (force)