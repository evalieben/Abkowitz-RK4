from ctypes import *
import matplotlib.pyplot as plt

def realtime_solve(dt,u,v,r,y,o,o_aim,X,Y,theta):
    ll = cdll.LoadLibrary   
    lib = ll("./Abkowitz.dll")
    data = []
    type_p_double = POINTER(c_double)
    lib.realtime_solve.argtypes = [c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double]	#设置输入参数类型
    lib.realtime_solve.restype  = type_p_double		#设置返回值类型
    f = lib.realtime_solve(dt,u,v,r,y,o,o_aim,X,Y,theta)
    for i in range(9):
        data.append(f[i])
    return data

def DtoR(D):
    ll = cdll.LoadLibrary   
    lib = ll("./Abkowitz.dll")
    lib.DegToRad.argtypes = [c_double]	#设置输入参数类型
    lib.DegToRad.restype  = c_double		#设置返回值类型
    return lib.DegToRad(D)

time = 0
dt = 0.1
rudder_angle = 15 #deg
data_t = [1,0,0,0,0,DtoR(rudder_angle),0,0,0]
t0 = 0
t1 = 200
X = []
Y = []
for t in range(int(t0/dt),int(t1/dt)):
    time = dt*t
    data_t = realtime_solve(dt,data_t[0],data_t[1],data_t[2],data_t[3],data_t[4],data_t[5],data_t[6],data_t[7],data_t[8])
    X.append(data_t[6])
    Y.append(data_t[7])
    #print(time,data_t)
    plt.ion()
    plt.scatter(X, Y,c = 'red')
    ax = plt.gca()
    ax.set_aspect(1)
    plt.pause(0.01)
    plt.clf()
    plt.show()
