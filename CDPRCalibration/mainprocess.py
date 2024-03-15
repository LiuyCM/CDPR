#Created by Liu yifan 2023-06-30
#Connect me by email: liu.y.cm@m.titech.ac.jp
#Last modified time: 2023-09-21

#import numpy as np
import matplotlib.pyplot as plt
import time_utils as tu
import CDPR_fuction as CF
import parameters_of_simplified_CDPR_model as PoSCDPRM
import none_linear_least_square as NLLS
#import geometric_elements as ge
import simulation_utils as su
import move_platform as MPlat
import file_io as fio
import error_effected as EE
#import post_process as pp
#import draw_utils as du
#from scipy.optimize import minimize

class MainProcess:
#Create a CDPR target
    #-----------------Function Start------------------
    def __init__(self,SAMNUM=None) -> None:
        if SAMNUM == None or not isinstance(SAMNUM,int):
            self.SAMNUM = 36
        else:
            self.SAMNUM = SAMNUM
        self.fig = plt.figure()
        time_marker = tu.time_stamp()
        #PreProcess
        self.MP = MPlat.move_platform()
        self.MP.initial(width=86,height=49,length = 287)
        #第一个CDPR用来保存真实的CDPR值，第二个保存设计的值
        self.CDPR = PoSCDPRM.CDPR(self.MP)
        self.CDPR1 = PoSCDPRM.CDPR(self.MP)
        self.EA = EE.errorAna(self.MP,self.CDPR,self.CDPR1)
        print("PreProcess takes " + tu.time_diff_str(time_marker))
        self.Create_plot()
        self.fio = fio.file()
    #------------------Function End-------------------

    def CDPR_initial(self):
        time_marker = tu.time_stamp()
        self.CDPR.initial(error = True)
        self.CDPR1.initial(error = False)
        #check the geometric errors
        #CDPR.debug_check_position()
        print("CDPR initialization takes " + tu.time_diff_str(time_marker))

#Create a 3D plot and axes
    #-----------------Function Start------------------
    def Create_plot(self):
        self.ax0 = self.fig.add_subplot(111,projection = '3d')
        self.ax_group = [self.ax0]
    #------------------Function End-------------------

#Simulate the calibration process
    #-----------------Function Start------------------
    def calibration(self,calibration_times = None):
        print("starting calibration process")
        time_marker = tu.time_stamp()
        #Do the calibration
        if calibration_times == None:
            count = 1
        else:
            count = calibration_times
        self.postProcesses = []
        while(count>0):
            postProcess = NLLS.NLLS_calibration(self.CDPR,ax = self.ax0,sam_num = self.SAMNUM,do_post_process = 1)
            self.postProcesses.append(postProcess)
            self.CDPR.after_calibration(postProcess)
            count -= 1
        print("Calibration takes time is " + tu.time_diff_str(time_marker))
        #------------------Function End-------------------


#Show the result of the calibration
    def show_result(self):
        #-----------------Function Start------------------
        #Do the post process

        for postProcess in self.postProcesses:
            time_marker = tu.time_stamp()
            postProcess.show_result()
            postProcess.print_accuracy()
            print("Post process takes time is " + tu.time_diff_str(time_marker))
            saveOrnot = input("Do you want to save the result? (y/n)")
            if saveOrnot == "y":
                postProcess.save_to_file()
        #------------------Function End-------------------


#Draw the initial CDPR position
    def draw(self):
        #-----------------Function Start------------------
        #self.CDPR.draw(self.ax0,mode = "complex")
        #CDPR.draw_workspace(ax1)
        #Set the ticks of the axes
        #Beware that the setting of the ticks should be after the drawing
        self.ax0.set_xlabel('X')
        self.ax0.set_ylabel('Y')
        self.ax0.set_zlabel('Z')
        self.ax0.set_box_aspect([1, 1, 1]) 
        self.ax0.set_xticks(self.CDPR.ticks.get_xticks())
        self.ax0.set_yticks(self.CDPR.ticks.get_yticks())
        self.ax0.set_zticks(self.CDPR.ticks.get_zticks())
        self.ax0.set_xticklabels(self.CDPR.ticks.get_xticks())
        self.ax0.set_yticklabels(self.CDPR.ticks.get_yticks())
        self.ax0.set_zticklabels(self.CDPR.ticks.get_zticks())
        plt.show()
        #------------------Function End-------------------


#visualize calibration result
#-----------------Function Start------------------
    def visualize_calibration_result(self,quick = None,is2D = None):
        if quick == None:
            quick = False
        else:
            quick = True
        if is2D == None:
            is2D = False
        else:
            is2D = True
        time_marker = tu.time_stamp()
        fig = plt.figure()
        ax0 = fig.add_subplot(221,projection = '3d')
        ax1 = fig.add_subplot(222,projection = '3d')
        ax2 = fig.add_subplot(223,projection = '3d')
        ax3 = fig.add_subplot(224,projection = '3d')
        su.error_evaluation(self.CDPR,fig,ax0,mode = "True",quick=quick,is2D=is2D)
        su.error_evaluation(self.CDPR,fig,ax1,mode = "Before",quick=quick,is2D=is2D)
        su.error_evaluation(self.CDPR,fig,ax2,mode = "After",quick=quick,is2D=is2D)
        su.error_evaluation(self.CDPR,fig,ax3,mode = "After2D",quick=quick,is2D=is2D)
        print("Visualizing calibration result takes " + tu.time_diff_str(time_marker , second = True))
        plt.show()
#------------------Function End-------------------


#模拟输入一个控制坐标后进行采样的过程
#输入：控制坐标（一个六维数组，分别表示xyz坐标和绕这三个轴旋转的角度）
#输出：各个绳索的受力和绳索的长度
#-----------------Function Start------------------
    def sample_simulation(self,pose):
        if len(pose) != 6:
            if len(pose) == 3:
                pose += [0,0,0]
            else:
                raise ValueError("The length of the pose should be 3 or 6")
        self.MP.set_pose(pose)
        self.MP.update_base_points()
        #time_markert = tu.time_stamp()
        B_fin,Ai = CF.IK(self.CDPR,self.MP,pulley=True,fin=True)
        #print("IK time:"+tu.time_diff_str(time_markert))
        #input("pause")
        self.MP.set_fin_points(B_fin)
        self.CDPR.set_cable_attach_points(Ai)
        CF.WS(self.CDPR,self.MP)

#-----------------Function End------------------
    def update_CDPR(self):
        self.CDPR.update_pulleys(self.MP)
        self.CDPR.update_force_sensors(self.MP)
            


#根据输入模式进行不同的操作
#这里是老代码，已经被重构了，保留仅为留档
'''
#-----------------Function Start------------------
    def run(self,mode = None):
        if mode == None:
            mode = "calibration"
        if mode == "calibration":
            self.CDPR_initial()
            self.calibration()
            self.show_result()
            self.draw_initial()
            ans = input("Do you want to visualize the calibration result? (y/n)")
            if ans == "y":
                is2 = input("Do you want to visualize the calibration result in 2D? (2/3/all)")
                if is2 == "2":
                    self.visualize_calibration_result(is2D=True)
                elif is2 == "3":
                    self.visualize_calibration_result()
                elif is2 == "all":
                    self.visualize_calibration_result()
                    self.visualize_calibration_result(is2D=True)
                else:
                    pass
            elif ans == "debug":
                self.visualize_calibration_result(quick=True)
                self.visualize_calibration_result(quick=True,is2D=True)
            else:
                print("Calibration is Over")
        elif mode == "debug":
            pass
        else:
            print("Wrong mode!")

#------------------Function End-------------------
'''

