#Created by Liu yifan 2023-09-28
#Connect me by email: liu.y.cm@m.titech.ac.jp
#Last modified time: 2024-02-08


import numpy as np
import draw_utils as du
import geometric_elements as ge
from scipy.optimize import least_squares
from scipy.optimize import approx_fprime
import simulation_utils as su
import post_process as pp
import time_utils as tu
import CDPR_fuction as CF
import jacobian as jac
import identification_matrix as im

class sample_info():
    def __init__(self,pose,lengths,forces):
        self.pose = pose
        self.length = lengths
        self.force = forces

    def check(self):
        print("pose:",self.pose)
        for i in range(1,5):
            print("length",i,":",self.length[i])
            print("force",i,":",self.force[i])
        

class NLLS():
    def __init__(self,CDPR,CDPR1,MP,samnum,sampledata = None):
        self.CDPR = CDPR
        self.CDPR1 = CDPR1
        self.MP = MP
        self.samnum = samnum
        #这里用的是基于条件数优化后的采样点默认16个点
        self.samnumsqrt = 4
        self.sam_info = []
        self.sam_true_points = []
        self.known = []
        self.unknown = []
        self.bound_lower = []
        self.bound_upper = []
        self.unknown_true = []
        self.n =  samnum * 6 + 12 + 4 #the number of unknowns
        self.m = samnum * 8 #the number of knowns
        self.result = None

        if sampledata == None:
            self.mp_inital=self.CDPR.move_platfrom_initial_position
            self.move_to([self.mp_inital.x,self.mp_inital.y,self.mp_inital.z])
            self.CDPR.initial_encoders()
            self.sam_init_points = self.generate_sampling()
            #self.sam_init_points = im.debug_sampling()
            self.sam_init_poses = []
            for i in range(len(self.sam_init_points)):
                self.sam_init_poses.append(self.sam_init_points[i]+[0,0,0])
        else:
            self.sampledata = sampledata
            self.move_to(sampledata.init_position)
            #忽视了弹性
            self.CDPR.initial_encoders(sampledata.init_length)
            self.sam_init_points = sampledata.point_position
            self.sam_init_poses = []
            for i in range(len(self.sam_init_points)):
                self.sam_init_poses.append(self.sam_init_points[i]+[0,0,0])
                info = sample_info(self.sam_init_poses[i],sampledata.read_relative_lengths(i),sampledata.cable_tension[i])
                self.sam_info.append(info)
            self.calibration_initial()
            self.debug_sample()
            input("pause")

        
    def move_to(self,pose,isRela = None,CDPR = None):
        if isRela is None or isRela != True:
            if len(pose) != 6:
                if len(pose) == 3:
                    pose += [0,0,0]
                else:
                    raise ValueError("The length of the pose should be 3 or 6")
        else:
            if len(pose) != 6:
                if len(pose) == 3:
                    pose += [0,0,0]
                    pose[0]+=self.mp_inital.x
                    pose[1]+=self.mp_inital.y
                    pose[2]+=self.mp_inital.z
                else:
                    raise ValueError("The length of the pose should be 3 or 6")
        self.MP.set_pose(pose)
        if CDPR is None:
            CF.WS(self.CDPR,self.MP)
            self.CDPR.update_pulleys(self.MP)
            self.CDPR.update_force_sensors(self.MP)
        else:
            CF.WS(CDPR,self.MP)
            CDPR.update_pulleys(self.MP)
            CDPR.update_force_sensors(self.MP)


    def generate_sampling(self):
        self.samx = np.array([0])
        self.samy = np.linspace(self.mp_inital.y-500,self.mp_inital.y+500,self.samnumsqrt)
        self.samz = np.linspace(self.mp_inital.z-1400,self.mp_inital.z+1400,self.samnumsqrt)
        self.samx,self.samy,self.samz = np.meshgrid(self.samx,self.samy,self.samz)
        self.samx = self.samx.reshape(-1)
        self.samy = self.samy.reshape(-1)
        self.samz = self.samz.reshape(-1)
        self.sam = np.array([self.samx,self.samy,self.samz])
        self.sam = self.sam.T
        return self.sam.tolist()
    
    def sample_all(self,withError = None):
        poses = []
        if withError is not None and withError == True:
            for i in range(len(self.sam_init_poses)):
                poses.append(su.generate_point_random_error(self.sam_init_poses[i][:3]+[0,0,0]))
        else:
            poses = self.sam_init_poses
        for i in range(len(self.sam_init_poses)):
            info = self.simu_MP_at_pose(poses[i])
            self.sam_info.append(info)
            self.sam_true_points.append(info.pose[:3])
            #self.CDPR.check_sensors()
            #self.sam_info[i].check()
    #这个模块用来检测假设已知滑轮位置，逆运动学计算结果是否能符合测量值
    def debug_sample(self):
        self.move_to(self.sampledata.init_position)
        print(self.sampledata.init_position)
        print(self.CDPR.read_force_sensors())
        templeLength = self.CDPR.read_encoders()
        for j in range(1,5):
            templeLength[j] += self.sampledata.init_length[j]
        print(templeLength)
        print(" ")
        for i in range(len(self.sam_init_poses)):
            self.move_to(self.sam_init_poses[i])
            print(self.sam_init_poses[i])
            print(self.CDPR.read_force_sensors())
            templeLength = self.CDPR.read_encoders()
            for j in range(1,5):
                templeLength[j] += self.sampledata.init_length[j]
            print(templeLength)
            print(" ")

    def simu_MP_at_pose(self,pose,isRela = None):
        self.move_to(pose,isRela)
        info = sample_info(self.MP.get_pose(),self.CDPR.read_encoders(),self.CDPR.read_force_sensors())
        return info

    def calibration_initial(self,withError = None):
        if withError is None:
            self.initial_knownlist()
            self.initial_unknownlist()
        else:
            self.initial_knownlist(withError)
            self.initial_unknownlist()

    def initial_knownlist(self,withError = None):
        self.known = []
        if withError is not None and withError == True:
            print("Warning:sensor error exist")
            for i in range(len(self.sam_info)):
                np.random.seed(2177)
                for j in range(1,5):
                    self.known.append(su.generate_length_sensor_error(self.sam_info[i].length[j]))
                for j in range(1,5):
                    self.known.append(su.generate_force_sensor_error(self.sam_info[i].force[j]))
        else:
            for i in range(len(self.sam_info)):
                for j in range(1,5):
                    self.known.append(self.sam_info[i].length[j])
                for j in range(1,5):
                    self.known.append(self.sam_info[i].force[j])
        #print("k=",str(len(self.known)/8))

    #未知数 = [a,L0,P]
    #P 代表采样点的位姿 len(P) = 6*k
    #a CDPR中pulley的基准坐标 len(a) = 3*4 = 12
    #L0 CDPR中cable的初始长度 len(L0) = 4
    def initial_unknownlist(self):
        self.unknown = []
        self.unknown_true = []
        self.bound_lower = []
        self.bound_upper = []
        for i in range(1,5):
            self.unknown.append(self.CDPR1.cable_base_points[i].x)
            self.unknown.append(self.CDPR1.cable_base_points[i].y)
            self.unknown.append(self.CDPR1.cable_base_points[i].z)
            self.unknown_true.append(self.CDPR.cable_base_points[i].x)
            self.unknown_true.append(self.CDPR.cable_base_points[i].y)
            self.unknown_true.append(self.CDPR.cable_base_points[i].z)

            self.bound_lower.append(self.CDPR1.cable_base_points[i].x-100)
            self.bound_lower.append(self.CDPR1.cable_base_points[i].y-100)
            self.bound_lower.append(self.CDPR1.cable_base_points[i].z-100)
            
            self.bound_upper.append(self.CDPR1.cable_base_points[i].x+100)
            self.bound_upper.append(self.CDPR1.cable_base_points[i].y+100)
            self.bound_upper.append(self.CDPR1.cable_base_points[i].z+100)
        for i in range(1,5):
            self.unknown.append(self.CDPR1.cable_length_designed)
            self.unknown_true.append(self.CDPR.encoders[i].init_distance)
            self.bound_lower.append(self.CDPR1.cable_length_designed-20)
            self.bound_upper.append(self.CDPR1.cable_length_designed+20)
        for i in range(len(self.sam_init_poses)):
            self.unknown += self.sam_init_poses[i]
            self.unknown_true += self.sam_info[i].pose
            self.bound_lower += [self.sam_init_poses[i][0]-100,self.sam_init_poses[i][1]-50,self.sam_init_poses[i][2]-50,self.sam_init_poses[i][3]-0.0001,self.sam_init_poses[i][4]-5,self.sam_init_poses[i][5]-5]
            self.bound_upper += [self.sam_init_poses[i][0]+100,self.sam_init_poses[i][1]+50,self.sam_init_poses[i][2]+50,self.sam_init_poses[i][3]+0.0001,self.sam_init_poses[i][4]+5,self.sam_init_poses[i][5]+5]
        print(len(self.unknown))

    def from_list_to_variables(self,x,issimplified = False):
        cable_base_points = {}
        for i in range(4):
            cable_base_points[i+1] = ge.Point3D(x[i*3],x[i*3+1],x[i*3+2])
        cable_init_length = {}
        for i in range(4):
            cable_init_length[i+1] = x[12+i]
        if issimplified == False:
            sampling_poses = []
            for i in range((len(x)-16)//6):
                sampling_poses.append([x[16+i*6],x[16+i*6+1],x[16+i*6+2],x[16+i*6+3],x[16+i*6+4],x[16+i*6+5]])
            return cable_base_points,cable_init_length,sampling_poses
        else:
            return cable_base_points,cable_init_length

    def from_variables_to_list(self,cable_base_points,cable_init_length,sampling_poses):
        list_x = []
        for i in range(4):
            list_x.append(cable_base_points[i+1][0])
            list_x.append(cable_base_points[i+1][1])
            list_x.append(cable_base_points[i+1][2])
        for i in range(4):
            list_x.append(cable_init_length[i+1])
        for i in range(len(sampling_poses)):
            list_x += sampling_poses[i]
        return list_x


    def _identification_matrix_norm(self,poses):
        print("calculating identification matrix...")
        time_markerIDM = tu.time_stamp()
        condition_numbers = []
        #jacs = []
        jax_y = -np.eye(8)
        jax_y_4 = -np.eye(4)
        column_pulley_x = [0,3,6,9]
        column_pulley_yz = [1,2,4,5,7,8,10,11]
        column_ini_cable = [12,13,14,15]
        column_EE_x = [16,20,21]
        column_EE_yz = [17,18,19]
        column_EE_po = [16,17,18]
        for i in range(len(poses)):
            jac_x = self._jac_diff(poses[i])
            needed_column = column_EE_x + column_EE_yz
            needed_column = sorted(needed_column)
            jac_x = jac_x[:,needed_column]
            #jac_x = jac_x[:4,:] #cable length
            jac_x = jac_x[-4:, :] #cable tension
            #jacs.append(jac_x)
            HP = np.linalg.inv(jax_y_4).dot(jac_x)
            idm = np.linalg.pinv(HP)
            #U, S, VT = np.linalg.svd(idm)
            #print("norm HP:",np.linalg.norm(HP))
            #print(HP)
            #print("norm idm:",np.linalg.norm(idm))
            #print(idm.T)
            condition_number = np.linalg.norm(HP)*np.linalg.norm(idm)
            print("for pose",i,"condition number is",condition_number)
            #print(S)
            #print(S[0]/S[-1])
            #input("pause")
            condition_numbers.append(condition_number)
        print("calculation finished!")
        print("sci_jac time:",tu.time_stamp()-time_markerIDM)
        return condition_numbers

    #use numerical method to calculate jacobi matrix
    def _jac_diff(self,pose = None):
        #print("calculating jacobi matrix...")
        #time_markerJAC = tu.time_stamp()
        #default epsilon = 1.49e-08
        if pose == None:
            sci_jac = approx_fprime(self.unknown,self.residucal_4_idm)
        else:
            cable_base_points,cable_init_length,_ = self.from_list_to_variables(self.unknown,issimplified = False)
            x_list = self.from_variables_to_list(cable_base_points,cable_init_length,[pose])
            sci_jac = approx_fprime(x_list,self.residucal_4_idm)
        #print("calculation finished!")
        #print("sci_jac time:",tu.time_stamp()-time_markerJAC)
        return sci_jac
    
    #use analytical method to calculate jacobi matrix
    def _jac_ana(self):
        print("the number of unknowns:",len(self.unknown))
        print("the number of knowns:",len(self.known))
        print("calculating jacobi matrix...")

        time_markerJAC = tu.time_stamp()
        cable_base_points,cable_init_length,sampling_poses = self.from_list_to_variables(self.unknown,issimplified = False)
        self.CDPR1.set_cable_base_points(cable_base_points)
        self.CDPR1.set_initial_length(cable_init_length)
        #for least square method, the jacobian matrix should be m*n, m is the number of variables, n is the number of equations
        #in this case, m = len(unknown), n = len(known)
        #First we initialize the jacobian matrix
        my_jac = np.zeros((len(self.known),len(self.unknown)))
        #the jacobian matrix is divided into two parts, one is the residuals of the cable length, the other is the residuals of cable tension
        #the length part is easy to calculate, the tension part is a little bit complicated
        jac_tension = np.zeros((len(self.known)//2,len(self.unknown)))
        jac_length = np.zeros((len(self.known)//2,len(self.unknown)))
        for i in range(len(sampling_poses)):
            jac_l_1,jac_l_2,jac_l_3 = jac.inv_jac_length(self.CDPR1,self.MP,sampling_poses[i])
            jac_length[i*4:i*4+4,0:12] = jac_l_1
            jac_length[i*4:i*4+4,12:16] = jac_l_2
            jac_length[i*4:i*4+4,i*6+16:i*6+22] = jac_l_3
        print("calculation finished!")
        print("sci_jac time:",tu.time_stamp()-time_markerJAC)
        #it should return my_jac, but for now, we only return jac_length
        return jac_length

    def residual(self,params):
        cable_base_points,cable_init_length,sampling_poses = self.from_list_to_variables(params,issimplified = False)
        self.CDPR1.set_cable_base_points(cable_base_points)
        self.CDPR1.set_initial_length(cable_init_length)
        #self.CDPR1.debug_check_position()
        resi = np.array([])
        cof_N = 1
        cof_NM = 1e-3
        cof_cons = 5
        for i in range(len(sampling_poses)):
            lengths,forces,balance_resi = CF.IC(self.CDPR1,self.MP,sampling_poses[i])
            for j in range(4):
                resi = np.append(resi,lengths[j+1]-self.known[i*8+j])
            for j in range(4):
                #this is the coefficient of force residual
                resi = np.append(resi,(forces[j]-self.known[i*8+j+4])*cof_N)
            resi = np.append(resi,balance_resi[0]*cof_N*cof_cons)
            resi = np.append(resi,balance_resi[1]*cof_NM*cof_cons)
            resi = np.append(resi,balance_resi[2]*cof_NM*cof_cons)
            resi = np.append(resi,0)
        ans = 0.5*resi.dot(resi.T)
        print("total residual",ans)
        return resi

    #this is a simplified version of residual function, which is used for identification matrix
    def residucal_4_idm(self,params):
        cable_base_points,cable_init_length,sampling_poses = self.from_list_to_variables(params,issimplified = False)
        self.CDPR1.set_cable_base_points(cable_base_points)
        self.CDPR1.set_initial_length(cable_init_length)
        #self.CDPR1.debug_check_position()
        resi = np.array([])
        cof_F = 10 #dN(real:140mN)
        for i in range(len(sampling_poses)):
            lengths,forces,_ = CF.IC(self.CDPR1,self.MP,sampling_poses[i])
            for j in range(4):
                resi = np.append(resi,lengths[j+1]-self.known[i*8+j])
            for j in range(4):
                #this is the coefficient of force residual
                resi = np.append(resi,(forces[j]-self.known[i*8+j+4])*cof_F)
        return resi
    

    def calibration_main(self,updateCDPR = None):
        time_markerCA = tu.time_stamp()
        result = least_squares(self.residual,self.unknown,method="trf",verbose=2,bounds = (self.bound_lower,self.bound_upper))
        #result = least_squares(self.residual,self.unknown,method="lm",verbose=2)
        
        print("Calibration takes "+tu.time_diff_str(time_markerCA))
        print(result)
        print(self.residual(self.unknown))
        print(self.residual(result.x))
        print(self.residual(self.unknown_true))
        print("-------------------------Result---------------------------------")

        cable_base_points,cable_init_length,sampling_points = self.from_list_to_variables(result.x,issimplified = False)
        self.result = (cable_base_points,cable_init_length,sampling_points)
        if updateCDPR is not None:
            updateCDPR.set_cable_base_points(cable_base_points)
            updateCDPR.initial_pulleys()
            updateCDPR.set_initial_length(cable_init_length)
        for i in range(1,5):
            print("Pulley: "+str(i)+" = ["+str(cable_base_points[i].x)+","+str(cable_base_points[i].y)+","+str(cable_base_points[i].z)+"]")
        for i in range(1,5):
            print("Cable: "+str(i)+" = "+str(cable_init_length[i]))
        for i in range(len(sampling_points)):
            print("MP: "+str(i)+" = "+str(sampling_points[i]))
        pulleyres =[]
        for i in range(2,5):
            pulleyres.append(cable_base_points[i].y-cable_base_points[1].y)
            pulleyres.append(cable_base_points[i].z-cable_base_points[1].z)

        print("-------------------------True---------------------------------")
        cable_base_pointsT,cable_init_lengthT,sampling_pointsT = self.from_list_to_variables(self.unknown_true,issimplified = False)
        for i in range(1,5):
            print("Pulley: "+str(i)+" = ["+str(cable_base_pointsT[i].x)+","+str(cable_base_pointsT[i].y)+","+str(cable_base_pointsT[i].z)+"]")
        for i in range(1,5):
            print("Cable: "+str(i)+" = "+str(cable_init_lengthT[i]))
        for i in range(len(sampling_pointsT)):
            print("MP: "+str(i)+" = "+str(sampling_pointsT[i]))
        pulleyresT =[]
        for i in range(2,5):
            pulleyresT.append(cable_base_pointsT[i].y-cable_base_pointsT[1].y)
            pulleyresT.append(cable_base_pointsT[i].z-cable_base_pointsT[1].z)
        print("--------------------------Compare--------------------------------")
        pulleybefor =[]
        for i in range(1,5):
            if i > 1:
                pulleybefor.append(self.CDPR.cable_base_points_designed[i].y-self.CDPR.cable_base_points_designed[1].y)
                pulleybefor.append(self.CDPR.cable_base_points_designed[i].z-self.CDPR.cable_base_points_designed[1].z)
            print("Pulley: "+str(i)+" = ["+str(self.CDPR.cable_base_points_designed[i].x-self.CDPR.cable_base_points_designed[1].x)+","+str(self.CDPR.cable_base_points_designed[i].y-self.CDPR.cable_base_points_designed[1].y)+","+str(self.CDPR.cable_base_points_designed[i].z-self.CDPR.cable_base_points_designed[1].z)+"]")
        for i in range(1,5):
            print("Pulley: "+str(i)+" = ["+str(cable_base_points[i].x-cable_base_points[1].x)+","+str(cable_base_points[i].y-cable_base_points[1].y)+","+str(cable_base_points[i].z-cable_base_points[1].z)+"]")
        for i in range(1,5):
            print("Pulley: "+str(i)+" = ["+str(cable_base_pointsT[i].x-cable_base_pointsT[1].x)+","+str(cable_base_pointsT[i].y-cable_base_pointsT[1].y)+","+str(cable_base_pointsT[i].z-cable_base_pointsT[1].z)+"]")
        acc = 0
        err = 0
        for i in range(len(pulleyres)):
            acc += abs(pulleyresT[i]-pulleyres[i])
            err += abs(pulleyresT[i]-pulleybefor[i])
        print("accuracy before calibration: "+str(err/6))
        print("accuracy after calibration: "+str(acc/6))
        print("error change: "+str(round(abs(acc/err)*100,2))+"%")












