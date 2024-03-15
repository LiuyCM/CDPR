#Created by Liu yifan 2023-11-12
#Connect me by email: liu.y.cm@titech.ac.jp
#Last modified time: 2023-11-12
#This file help to output some data
import time
import numpy as np
class file():
    def __init__(self,filename = None):
        self.starttime = time.time()
        if filename is None:
            #filename should follow the format "month-day-year-hour-minute-second.txt"
            #and the file will be saved in a folder named "record"
            #example: 11-12-2023:12-30-00.txt
            self.filename = time.strftime("record/%m-%d-%Y_%H-%M-%S") + ".txt"
        else:
            self.filename = "record/" + filename + ".txt"
        #create a file
        self.file = open(self.filename,"w")
        self.file.close()

    def save_array_to_file(self,array,file_name = None):
        if file_name is None:
            file_name = self.filename
        with open(file_name,"a") as file:
            np.savetxt(file,array)

    def save_jac_to_file(self,jac_array,file_name = None,wid=None,pre=None,lenres = None):
        if file_name is None:
            file_name = self.filename
        if wid is None:
            wid = 16
        if pre is None:
            pre = 9
        if lenres is None:
            lenres = 8
        with open(file_name, 'a') as f:
            #first, output some information about the jacobian matrix like the size of the matrix
            f.write("Jacobian matrix size: " + str(jac_array.shape[0]) + " * " + str(jac_array.shape[1]) + '\n')
            line = ' '.join('{: 0{width}d}'.format(i+1,width=wid) for i in range(jac_array.shape[1]))
            f.write(line + '\n')
            f.write('-------------------------------------------------------------------------------------------\n')
            counter_t = 0
            for row in jac_array:
                counter_t += 1
                #To control the precision of the output, we use the format .5e
                line = ' '.join('{: >{width}.{precision}e}'.format(x,width=wid,precision=pre) for x in row)
                f.write(line + '\n')
                if counter_t % lenres == 0:
                    f.write('\n')
                    counter_t = 0

    def save_points_to_file(self,points,file_name = None):
        if file_name is None:
            file_name = self.filename
        with open(file_name,"a") as file:
            for point in points:
                file.write(str(point[0]) + " " + str(point[1]) + " " + str(point[2]) + '\n')
            file.write('\n')



