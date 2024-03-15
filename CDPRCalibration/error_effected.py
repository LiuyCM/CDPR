#Created by Liu yifan 2024-01-03
#Connect me by email: liu.y.cm@m.titech.ac.jp
#Last modified time: 2024-01-15
#Õâ¸öº¯ÊýÓÃÀ´¼Eé¼¸ºÎ²ÎÊýÎó²ûÒÔÓÚ×ûòÕÄ©¶ËÖ´ÐÐÆ÷Î»ÖÃµÄÓ°ÏE
import numpy as np
import geometric_elements as ge
import CDPR_fuction as CF
import time_utils as tu
from scipy.optimize import minimize
from scipy.optimize import least_squares

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

class errorAna():
    def __init__(self,Moveplatform,CDPR,CDPR1):
        self.MP = Moveplatform
        self.CDPRWithError = CDPR
        self.CDPRWithoutError = CDPR1
        self.CDPRWithoutError.initial(error = False)
        self.stiffunit = 59800 #N/m
        self.LPulToPul = {} #mm
        self.LPulToPul[1] = 6309
        self.LPulToPul[2] = 6503
        self.LPulToPul[3] = 6503
        self.LPulToPul[4] = 6309
        self.ld = {}
        self.Td = {}
        self.Ld = {}
        self.pdf = {}
        self.g_MP = np.array([0,0,-self.MP.get_mass()*9.8])
    
     #¸Ãº¯ÊýÓÃÀ´¸EÂº¬Îó²ûÑÄCDPRÄ£ÐÍ£¬½ÓÊÜÊäÈE¬ÊäÈEª»¬ÂÖÏà¶ÔÉè¼ÆÎ»ÖÃµÄÎó²û¿¬×ÖµäÐÎÊ?
    def set_CDPR(self,errors = None):
        if errors is None:
            errors = {}
            errors[1] = [0,0,0]
            errors[2] = [50,0,0]
            errors[3] = [50,0,0]
            errors[4] = [0,0,0]
        self.CDPRWithError.initial(error = True,set_errors = errors)
        
    #¼Eéº¬Îó²ûÑÄCDPRµÄÕæÖµ
    def get_pulleys_error(self):
        for i in range(1,5):
            print(self.CDPRWithError.cable_base_points[i])

    #·â×°Õû¸öÁ÷³Ì,²»°E¬ÉèÖÃÄ£ÐÍÎó²E
    def pipeline(self,pose):
        pose = self.posecheck(pose)
        self.controlled_LandT(pose)
        self.update_Ld()
        poseWithoutATFB,cable_forces,flag = self.search_pose(pose)
        #input("finish poseWithoutATFB")
        #print(poseWithoutATFB)
        #print(cable_forces)
        #print(self.average_tension(cable_forces))
        if flag == -1:
            print("The pose is not reachable")
            return poseWithoutATFB,poseWithoutATFB
        poseWithATFB = self.ATFB_pose(poseWithoutATFB,cable_forces)
        return poseWithoutATFB, poseWithATFB

    #Step1£¬Calculate the controlled cable lengths and cable tension under designed parameters
    def controlled_LandT(self,pose):
        self.move_to(pose,self.CDPRWithoutError)
        force = self.CDPRWithoutError.read_force_sensors()
        length = self.CDPRWithoutError.read_encoders()
        print("The controlled force is :",force)
        print("The controlled length is :",length)
        self.ld = length
        self.Td = force
        B_fin,Ai = CF.IK(self.CDPRWithoutError,self.MP,pulley=True,fin=True)
        u = self.unit_vectors(Ai,B_fin)
        b = {}
        omega = {}
        for i in range(1,5):
            u[i] = np.array([u[i][1],u[i][2]])
            b[i] = np.array([B_fin[i][1]-pose[1],B_fin[i][2]-pose[2]])
            omega[i] = np.array([u[i][0],u[i][1],np.cross(b[i],u[i])])
        pdf_list = np.vstack((omega[1],omega[2],omega[3]))
        pdf_list = -1*np.linalg.inv(pdf_list.T)
        pdf_list = pdf_list.dot(omega[4].T)
        pdf_list = pdf_list.T.tolist() + [1]
        for i in range(1,5):
            self.pdf[i] = pdf_list[i-1]
        print(self.pdf)

    #Step2:Take the cable as spring and calculate the rest length
    def update_Ld(self):
        for i in range(1,5):
            self.Ld[i] = self.ld[i] + self.LPulToPul[i] - self.Td[i] * ((self.ld[i] + self.LPulToPul[i])/self.stiffunit)
            print("The rest length of cable",i,"is:",self.Ld[i])

    #Step3: Search the equilibrium pose of the end-effector under geometric errors. 
    def search_pose(self,pose):
        time_markerSP = tu.time_stamp()
        sigma = 2e-1
        kf = 1e-1
        #Step3.1: find a orientaion of the end-effector make the moment of the force is zero
        orientation,ela_forces = self.best_orientataion(pose)
        #Step3.2: calculate the net force of the end-effector
        net_force = self.calcu_net_force(ela_forces)
        #print("netforce:",net_force)
        norm_net_force = np.linalg.norm(net_force)
        while norm_net_force > sigma:
            #print("start norm_net_force:",norm_net_force)
            newpose = pose[:3] + kf * net_force
            #newpose = newpose.tolist() + orientation.tolist()
            new_orientation,new_ela_forces = self.best_orientataion(newpose)
            new_net_force = self.calcu_net_force(new_ela_forces)
            norm_new_net_force = np.linalg.norm(new_net_force)
            
            #print("kf:",kf)
            #print("norm_new_net_force:",norm_new_net_force)
            
            if norm_new_net_force > norm_net_force:
                kf = kf*0.5
                #print("step too large,shrinke the kf",kf)
            else:
                pose = newpose
                orientation = new_orientation
                ela_forces = new_ela_forces
                net_force = new_net_force
                norm_net_force = norm_new_net_force
                kf = kf*1.1
                #print("step acceptable,update the pose",kf)
                print("norm_net_force:",norm_net_force)
            if kf < 1e-8:
                print("step too small")
                print("---------------Debug Info---------------------------")
                print("pose",pose)
                print("newpose",newpose)
                print("orientation:",orientation)
                print("new_orientation:",new_orientation)
                print("net_force:",net_force)
                print("new_net_force:",new_net_force)
                print("norm_net_force",norm_net_force)
                print("norm_new_net_force",norm_new_net_force)
                print("move vector:",kf*new_net_force)
                print("---------------Debug Info---------------------------")
                #input("pause")
                if norm_new_net_force < 2*sigma:
                    return [newpose[0],newpose[1],newpose[2],new_orientation[0],new_orientation[1],new_orientation[2]],new_ela_forces,1
                else:
                    return [newpose[0],newpose[1],newpose[2],new_orientation[0],new_orientation[1],new_orientation[2]],new_ela_forces,-1
        #print("Find actually pose takes " + tu.time_diff_str(time_markerSP))
        return [pose[0],pose[1],pose[2],orientation[0],orientation[1],orientation[2]],ela_forces,1

    # Step4: adjust the length of the cables according to the average tension
    def ATFB_pose(self, pose,cable_forces):
        sigma = 2e-1
        Kten = 1
        time_markerAP = tu.time_stamp()
        avg_tension = self.average_tension(cable_forces)
        tensionFB = avg_tension - 100
        #print("avg_tension:",avg_tension)
        #print("tensionFB",tensionFB)
        #print("Start Adjust")
        while abs(tensionFB) > sigma:
            #Å£¶Ù·¨update Ld
            oldLd = {}
            for i in range(1,5):
                oldLd[i] = self.Ld[i]
                self.Ld[i] = self.Ld[i] + (Kten * self.pdf[i] * tensionFB * ((self.ld[i] + self.LPulToPul[i])/self.stiffunit))
            #print("cable length before adjust:",oldLd)
            #print("cable length after adjust:",self.Ld)
            poseWithATFB,cable_forces,flag = self.search_pose(pose)
            if flag == -1:
                for i in range(1,5):
                    self.Ld[i] = oldLd[i]
                Kten = Kten * 0.5
            else:
                new_avg_tension = self.average_tension(cable_forces)
                new_tensionFB = new_avg_tension - 100
                print("new_avg_tension:",new_avg_tension)
                if abs(new_tensionFB) < abs(tensionFB):
                    #print("result acceptable")
                    tensionFB = new_tensionFB
                    avg_tension = new_avg_tension
                    pose = poseWithATFB
                    Kten = Kten * 1.1
                else:
                    #print("result unacceptable")
                    for i in range(1,5):
                        self.Ld[i] = oldLd[i]
                    Kten = Kten * 0.5
                #input("Press Enter to continue...")
            #print("kten",Kten)
            
        #print("Adjust the length of the cables takes " + tu.time_diff_str(time_markerAP))
        return pose


        

    def average_tension(self,cable_tension_vectors):
        avgT = 0
        for i in range(1,5):
            tension = np.linalg.norm(cable_tension_vectors[i])
            avgT += abs(tension)
        return avgT/4

    def calcu_net_force(self,forces):
        net_force = np.array([0.0,0.0,0.0])
        for i in range(1,5):
            net_force += forces[i]
        return net_force+self.g_MP

    def best_orientataion(self,pose):
        #time_markerBO = tu.time_stamp()
        force_vector = {}
        def moment_residual(orientation):
            self.MP.set_pose([pose[0],pose[1],pose[2],orientation[0],orientation[1],orientation[2]])
            originP = np.array([pose[0],pose[1],pose[2]])
            B_fin,Ai = CF.IK(self.CDPRWithError,self.MP,pulley=True,fin=True)
            unit_force = self.unit_vectors(Ai,B_fin)
            self.MP.set_fin_points(B_fin)
            self.CDPRWithError.set_cable_attach_points(Ai)
            self.CDPRWithError.update_pulleys(self.MP)
            templength = self.CDPRWithError.read_encoders()
            elasticforce = self.calcu_elastic_force(templength)
            #print("length:",templength)
            #print("force:",elasticforce)
            moment = 0
            for i in range(1,5):
                force_vector[i] = elasticforce[i] * unit_force[i]
                moment += np.cross(force_vector[i],B_fin[i]-originP)
            return moment
        orient = np.array([0,0,0])
        bnds = ([-90,-90,-90],[90,90,90])
        result = least_squares(moment_residual,orient,ftol=1e-6, xtol=1e-6, gtol=1e-6,bounds=bnds,verbose=False)
        #print("Find best orientataion takes " + tu.time_diff_str(time_markerBO))
        return result.x,force_vector

    def posecheck(self,pose,isRela = None):
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
        return pose

    def move_to(self,pose,CDPR):
        self.MP.set_pose(pose)
        CF.WS(CDPR,self.MP)
        CDPR.update_pulleys(self.MP)
        CDPR.update_force_sensors(self.MP)

    def calcu_elastic_force(self,slength):
        eforce = {}
        for i in range(1,5):
            dl = slength[i] + self.LPulToPul[i] - self.Ld[i]
            if dl < 0:
                eforce[i] = 0
            else:
                eforce[i] = dl * self.stiffunit / (slength[i] + self.LPulToPul[i])
        return eforce

    def unit_vectors(self,A,B):
        u = {}
        for i in range(1,5):
            u[i] = (A[i] - B[i])/np.linalg.norm(A[i] - B[i])
        return u

