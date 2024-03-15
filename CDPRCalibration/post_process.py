#Created by Liu yifan 2023-05-27
#Connect me by email: liu.y.cm@titech.ac.jp
#Last modified time: 2023-06-07
#This file contains some functions that can be used to pocess the calibration result
#import numpy as np
import none_linear_least_square as NLLS
import datetime

#print the result beautifully
def pprintx(x,samp = None):
    print("cable base point:")
    for i in range(4):
        print("A"+str(i)+": "+str(x[i*3:i*3 + 3]))
    print("cable length: "+str(x[12:16]))
    num = len(x)
    print("sampling times:"+str((num-16)//3))
    if samp == None:
        print("sampling points:")
        for i in range((num-16)//3):
            print("P"+str(i+1)+": "+str(x[16+i*3:16+i*3+3]))

def fpprintx(x, samp=None):
    output = ""
    output += "cable base point:\n"
    for i in range(4):
        output += "A"+str(i)+": "+str(x[i*3:i*3 + 3]) + "\n"
    output += "cable length: "+str(x[12:16]) + "\n"
    num = len(x)
    output += "sampling times:"+str((num-16)//3) + "\n"
    if samp == None:
        output += "sampling points:\n"
        for i in range((num-16)//3):
            output += "P"+str(i+1)+": "+str(x[16+i*3:16+i*3+3]) + "\n"
    return output

def pprinty(y):
    print("encoder value:")
    for index in y.keys():
        print(str(index) + " : " +str(y[index]))

def value_resi(result):
    ans = 0
    for v in result:
        ans += v**2
    return ans

def generate_combinations(lst):
    if len(lst) == 0:
        return [[]]
    else:
        sub_combinations = generate_combinations(lst[1:])
        new_combinations = []
        for combination in sub_combinations:
            new_combinations.append([+lst[0]] + combination)
            new_combinations.append([-lst[0]] + combination)
        return new_combinations
#Useless
'''
def inverse_check(x,**kwargs):
    newx = x[:12].tolist()
    other_part_x = x[12:].tolist()
    can_inverse_x = newx[:5] + [newx[6]] + [newx[9]] + [newx[11]]
    print(can_inverse_x)
    inverse_resulte = generate_combinations(can_inverse_x)
    resi = NLLS.residual(x,**kwargs)
    min_v = value_resi(resi)
    print(str(min_v) + " is the initial value")
    best_nx = x
    for nx in inverse_resulte:
        nx = nx[:5] + [newx[5]] + [nx[5]] + newx[7:9] + [nx[6]] + [newx[10]] + [nx[7]] + other_part_x[:]
        resi = NLLS.residual(nx,**kwargs)
        vr = value_resi(resi)
        print(vr)
        if vr < min_v:
            min_v = vr(resi)
            best_nx = nx
    print(best_nx)
'''
#输入设计值或校准后值和真值，返回各分量的误差
def the_diff_between_two_params(params1,params2):
    if(len(params1) != len(params2)):
        return "error"
    else:
        diff = []
        for i in range(len(params1)):
            diff.append(params1[i] - params2[i])
        return diff


class post_processing:
    def __init__(self,result,x,y,x_true,timeused):
        self.result = result
        self.result_x = result.x
        self.x_init = x
        self.x_true = x_true
        self.error_init = the_diff_between_two_params(self.x_init,self.x_true)
        self.error_result = the_diff_between_two_params(self.result_x,self.x_true)
        self.y = y
        self.cable_base_points_init,self.cable_init_length = NLLS.from_list_to_vaiables(x,issimplified = True)
        self.cable_base_points_true,self.cable_init_length_true = NLLS.from_list_to_vaiables(x_true,issimplified = True)
        self.cable_base_points_result,self.cable_init_length_result = NLLS.from_list_to_vaiables(result.x,issimplified = True)
        self.cable_base_point2Ds_result = self.from_3DPoints_to_2DPoints()
        self.cable_base_points_error_init,self.cable_init_length_error_init = NLLS.from_list_to_vaiables(self.error_init,issimplified = True)
        self.cable_base_points_error_result,self.cable_init_length_error_result = NLLS.from_list_to_vaiables(self.error_result,issimplified = True) 
        self.timeused = timeused
        self.cable_base_points_accuracy, self.cable_length_accuracy = self.cal_accuracy()
        self.cable_base_points2D_accuracy, self.cable_length_accuracy = self.cal_accuracy(ignoreX=True)

    def from_3DPoints_to_2DPoints(self):
        cable_base_point2Ds = {}
        for index,point in self.cable_base_points_result.items():
            p = point.copy()
            p.x = 0
            cable_base_point2Ds[index] = p
        return cable_base_point2Ds
    def cal_result(self):
        for i in range(len(self.cable_base_points_init)):
            self.cable_base_points_error_init[i] = self.cable_base_points_true[i] - self.cable_base_points_init[i]
            self.cable_base_points_error_result[i] = self.cable_base_points_true[i] - self.cable_base_points_result[i]
        for i in range(len(self.cable_init_length)):
            self.cable_init_length_error_init[i] = self.cable_init_length_true[i] - self.cable_init_length[i]
            self.cable_init_length_error_result[i] = self.cable_init_length_true[i] - self.cable_init_length_result[i]

    def record_time(self,time):
        self.timeused = time

    def cal_accuracy(self,ignoreX = None):
        cable_base_points_accuracy = {}
        cable_base_points2D_accuracy = {}
        cable_length_accuracy = {}

        # 计算基准点的校准精度
        if ignoreX == None:
            for i in range(len(self.cable_base_points_init)):
                init_error = self.cable_base_points_error_init[i]
                result_error = self.cable_base_points_error_result[i]

                if init_error != 0:
                    accuracy = (result_error / init_error)
                else:
                    accuracy = float('inf')
                cable_base_points_accuracy[i] = accuracy
        else:
            for i in range(len(self.cable_base_points_init)):
                init_error = self.cable_base_points_error_init[i]
                result_error = self.cable_base_points_error_result[i]
                if init_error != 0 :
                    accuracy = result_error.div2D(init_error)
                else:
                    accuracy = float('inf')
                cable_base_points2D_accuracy[i] = accuracy

        # 计算长度的校准精度
        for i in range(len(self.cable_init_length)):
            init_error = self.cable_init_length_error_init[i]
            result_error = self.cable_init_length_error_result[i]

            if init_error != 0:
                accuracy = abs(result_error / init_error)
            else:
                accuracy = float('inf')

            cable_length_accuracy[i] = accuracy
        if ignoreX == None:
            return cable_base_points_accuracy, cable_length_accuracy
        else:
            return cable_base_points2D_accuracy, cable_length_accuracy


    def print_accuracy(self):
        print("-----------------------")
        print("the accuracy of the result is: ")
        for i in range(len(self.cable_base_points_init)):
            print("base point" + str(i) + "'s accuracy: " + "{:.3f}%".format(self.cable_base_points_accuracy[i]*100))
        for i in range(len(self.cable_base_points_init)):
            print("base point2D" + str(i) + "'s accuracy: " + "{:.3f}%".format(self.cable_base_points2D_accuracy[i]*100))
        for i in range(len(self.cable_init_length)):
            print("length" + str(i) + "'s accuracy: " + "{:.3f}%".format(self.cable_length_accuracy[i]*100))

    def show_result(self):
        print("-----------------------")
        print("the inital value is: ")
        pprintx(self.x_init,samp = 1)
        print("-----------------------")
        print("the true value is: ")
        pprintx(self.x_true,samp = 1)
        print("-----------------------")
        print("the result value is: ")
        pprintx(self.result_x,samp = 1)
        print("-----------------------")
        print("the result is: ")
        print(self.result)
        print("-----------------------")
        print("the error between init and true is: ")
        pprintx(self.error_init,samp = 1)
        print("-----------------------")
        print("the error between result and true is: ")
        pprintx(self.error_result,samp = 1)
        print("-----------------------")
        print("Calibration time is " + "{:.3f}".format(self.timeused) + " s")

    def save_to_file(self):
        filename = datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S.txt')
        with open(filename, 'w') as file:
            file.write("-----------------------\n")
            file.write("the accuracy of the result is: \n")
            for i in range(len(self.cable_base_points_init)):
                file.write("base point" + str(i) + "'s accuracy: " + "{:.3f}%".format(self.cable_base_points_accuracy[i]*100) + "\n")
            for i in range(len(self.cable_base_points_init)):
                file.write("base point2D" + str(i) + "'s accuracy: " + "{:.3f}%".format(self.cable_base_points2D_accuracy[i]*100) + "\n")
            for i in range(len(self.cable_init_length)):
                file.write("length" + str(i) + "'s accuracy: " + "{:.3f}%".format(self.cable_length_accuracy[i]*100) + "\n")
            
            file.write("-----------------------\n")
            file.write("the inital value is: \n")
            file.write(fpprintx(self.x_init,samp = 1) + "\n")
            file.write("-----------------------\n")
            file.write("the true value is: \n")
            file.write(fpprintx(self.x_true,samp = 1) + "\n")
            file.write("-----------------------\n")
            file.write("the result value is: \n")
            file.write(fpprintx(self.result_x,samp = 1) + "\n")
            file.write("-----------------------\n")
            file.write("the result is: \n")
            file.write(str(self.result) + "\n")
            file.write("-----------------------\n")
            file.write("the error between init and true is: \n")
            file.write(fpprintx(self.error_init,samp = 1) + "\n")
            file.write("-----------------------\n")
            file.write("the error between result and true is: \n")
            file.write(fpprintx(self.error_result,samp = 1) + "\n")
            file.write("-----------------------\n")
            file.write("Calibration time is " + str(self.timeused) + " s" + "\n")

