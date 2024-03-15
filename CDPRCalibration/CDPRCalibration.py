#Created by Liu yifan 2023-05-26
#Connect me by email: liu.y.cm@m.titech.ac.jp
#Last modified time: 2024-02-08
#开发版本对于原有函数和逻辑进行了大幅度修改，旧函数是没法正常工作的，如果要用老版本去以前semi里面找
import mainprocess as mp
import CDPR_fuction as CF
import numpy as np
import matplotlib.pyplot as plt
import time_utils as tu
#import parameters_of_simplified_CDPR_model as PoSCDPRM
#import none_linear_least_square as NLLS
#import geometric_elements as ge
import simulation_utils as su
#import post_process as pp
import draw_utils as du
import calibration_function as CalFun
import experiment_result as er

#from scipy.optimize import minimize

#记得最后把子类合并回去，不然就成屎山代码了
class debugMP(mp.MainProcess):
    def __init__(self,SAMNUM,sampledata = None):
        super().__init__(SAMNUM)
        self.simulateFlag = False
        if sampledata == None:
            self.simulateFlag = True
        else:
            self.sampledata = sampledata
            self.SAMNUM = sampledata.samnum

    def run(self,mode = None):
        if mode == None:
            mode = "simulation"
        if mode == "calibration":
            if self.simulateFlag == True:
                print("Please input the sample data")
                return 
            #self.sampledata.check()
            self.CDPR_initial()
            Calibration = CalFun.NLLS(self.CDPR,self.CDPR1,self.MP,self.SAMNUM,self.sampledata)
            Calibration.calibration_main(updateCDPR=self.CDPR1)
            #return 0
            #self.show_result()
            #du.draw_points(Calibration.sam_init_points,self.ax0,label = 'theoretical sampling poses')
            du.draw_points(Calibration.sam_true_points,self.ax0,color= 'c',label = 'sampling poses')
            Calibration.move_to(Calibration.sam_info[int(self.SAMNUM/2)].pose)
            du.draw_CDPR(CDPR=self.CDPR,MP=self.MP,ax=self.ax0,Pulley = True,Fin = True,showlegend=True)
            Calibration.move_to(Calibration.result[2][int(SAMNUM/2)],CDPR=self.CDPR1)
            #du.draw_CDPR(CDPR=self.CDPR1,MP=self.MP,ax=self.ax0,Pulley = True,Fin = True)
            self.draw()
            return 0
        elif mode == "simulation":
            #初始化随机的含误差CDPR模型
            #self.CDPR.debug_check_position()
            #以下为标准仿真程序
            self.CDPR_initial()
            Calibration = CalFun.NLLS(self.CDPR,self.CDPR1,self.MP,self.SAMNUM)
            Calibration.sample_all(withError=True)
            Calibration.calibration_initial(withError=True)
            #Calibration.calibration_initial()
            Calibration.calibration_main(updateCDPR=self.CDPR1)

            #du.draw_points(Calibration.sam_init_points,self.ax0,label = 'theoretical sampling poses')
            du.draw_points(Calibration.sam_true_points,self.ax0,color= 'c',label = 'sampling poses')
            Calibration.move_to(Calibration.sam_info[int(self.SAMNUM/2)].pose)
            du.draw_CDPR(CDPR=self.CDPR,MP=self.MP,ax=self.ax0,Pulley = True,Fin = True,showlegend=True)
            Calibration.move_to(Calibration.result[2][int(SAMNUM/2)],CDPR=self.CDPR1)
            #du.draw_CDPR(CDPR=self.CDPR1,MP=self.MP,ax=self.ax0,Pulley = True,Fin = True)
            self.draw()
            return 0
            
        
        elif mode == "error_analysis":
            #------------------------------------------------------
            #这个部分是估计滑轮误差对最终位置的影响。
            self.EA.set_CDPR()
            self.EA.get_pulleys_error()
            #points = [[0,300,500],[0,600,1500],[0,350,3000],[0,600,3000],[0,1200,3000],[0,1200,500]]
            points = [[0,400,3000]]
            for point in points:
                Pwithout, Pwith = self.EA.pipeline(point)
                #print(Pwithout)
                print(point)
                print(Pwith)
                Diserror = [0,0,0]
                totalserror = 0
                for i in range(3):
                    Diserror[i] = (Pwith[i]-point[i])
                RM=su.get_rotation_matrix(Pwith[3:])
                norm = np.array([-1,0,0])
                Rnorm = np.dot(RM,norm)
                yzRnorm = np.array([0,Rnorm[1],Rnorm[2]])
                yzRnorm = yzRnorm/np.linalg.norm(yzRnorm)
                print(Rnorm)
                Rangle = np.arccos(np.dot(Rnorm,norm))
                if abs(np.dot(Rnorm,norm)-1) < 1e-6:
                    Rangle = 0
                print(np.degrees(Rangle))
                orierror = (1000 * np.tan(Rangle)) * yzRnorm
                truetarget = Diserror + orierror
                print("Diserror is ",Diserror)
                print("orierror is ",orierror)
                print("truetarget is ",truetarget)
                totalserror = np.linalg.norm(truetarget)
                print("totalserror is ",totalserror)
            du.draw_points(points,self.ax0,color= 'c',label = 'controlled point')
            #du.draw_points([Pwithout[:3]],self.ax0,color= 'y',label = 'point without ATFB')
            du.draw_points([Pwith[:3]],self.ax0,color= 'b',label = 'point with ATFB')
            du.draw_line(points[-1], Pwithout[:3], self.ax0, color='y')
            du.draw_line(Pwithout[:3], Pwith[:3], self.ax0, color='r')
            du.draw_CDPR(CDPR=self.CDPR,MP=self.MP,ax=self.ax0,Pulley = True,Fin = True,showlegend=True)
            self.draw()
            return 0
            #------------------------------------------------------
        elif mode == "identification":
            #------------------------------------------------------
            #这个部分是计算条件数的，也就是在哪里采样效率高。
            pnum = 8
            samx = np.array([0])
            samy = np.linspace(Calibration.mp_inital.y-500,Calibration.mp_inital.y+500,pnum)
            samz = np.linspace(Calibration.mp_inital.z-1400,Calibration.mp_inital.z+1400,int(pnum*2))
            samx,samy,samz = np.meshgrid(samx,samy,samz)
            samx = samx.reshape(-1)
            samy = samy.reshape(-1)
            samz = samz.reshape(-1)
            sampoints = np.array([samx,samy,samz])
            sampoints = sampoints.T.tolist()
            samposes= []
            for i in range(len(sampoints)):
                samposes.append(sampoints[i]+[0,0,0])
            print("we need calculate ", 2*pnum**2,"points, pls wait patiently")
            test = Calibration._identification_matrix_norm(poses=samposes)
            
            # 将列表元素及其索引组合成元组列表
            indexed_lst = list(enumerate(test))
            # 根据元素的值对列表进行排序
            sorted_lst = sorted(indexed_lst, key=lambda x: x[1])
            # 取出最小的n个值的索引
            smallest_indices = [index for index, value in sorted_lst[:int(0.3*pnum**2)]]
            smallest_points = [sampoints[i] for i in smallest_indices]
            # 取出最大的n个值的索引
            sorted_lst = sorted(indexed_lst, key=lambda x: x[1], reverse=True)
            largest_indices = [index for index, value in sorted_lst[:int(0.3*pnum**2)]]
            largest_points = [sampoints[i] for i in largest_indices]
            
            #self.fio.save_jac_to_file(idmtest,lenres = 8)
            self.fio.save_points_to_file(smallest_points)
            #con = input("test has been saved. Continue?")
            #jactest = Calibration._jac_ana()
            du.draw_point_cloud(points=sampoints,values=test,ax=self.ax0,fig=self.fig)
            #du.draw_points(smallest_points,self.ax0,color= 'b',label = 'small sampling poses')
            #du.draw_points(largest_points,self.ax0,color= 'r',label = 'large sampling poses')
            '''
            plt.clf()
            self.fig = plt.figure()
            self.ax0 = self.fig.add_subplot(111,projection = '3d')
            du.draw_points(smallest_points,self.ax0,color= 'b',label = 'best sampling poses')
            self.draw()
            '''
        else:
            print("Wrong mode!")

experimentData = er.ExResult()
SAMNUM = 9
main = debugMP(SAMNUM,experimentData)
#main.run("simulation")
main.run("calibration")

