from ctypes import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import threading
import keyboard,time

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

def RtoD(R):
    return R*180/3.1415926

rudder_angle = 0 
o_aim = 0
stop = False
pause = False
def actuateRudder():
    global o_aim
    global stop
    global pause
    global free
    global zig
    global turn
    while(not stop):
        if(free == True):
            if keyboard.is_pressed('a'):
                    if(o_aim>=25):
                        pass
                    else:
                        o_aim+=0.5
            if keyboard.is_pressed('d'):
                if(o_aim<=-25):
                    pass
                else:
                    o_aim-=0.5
        if keyboard.is_pressed('q'):
            stop = True
        if keyboard.is_pressed('z'):
            zig = True
            turn = False
            free = False
            time.sleep(1)
        if keyboard.is_pressed('x'):
            turn = True
            zig = False
            free = False
        if keyboard.is_pressed('c'):
            free = True
            zig = False
            turn = False
            time.sleep(1)
        if keyboard.is_pressed('SPACE'):
            pause = bool(1-pause)
            time.sleep(1)
        time.sleep(0.01)

def zigzag(o_t,theta_t):
    global o_aim
    if theta_t-2*3.1415926*int(theta_t/(2*3.1415926)) > DtoR(o_aim):
        o_aim = -15
    else:
        o_aim = 15
def turn10():
    global o_aim
    o_aim = 10

Time = 0
dt = 0.1
data_t = [1,0,0,0,0,DtoR(rudder_angle),0,0,0]
X = []
Y = []
U = []
V = []
R = []
O = []
O_AIM = []
THETA = []
TIME = []
t = 0
zig = False
turn = False
free = True

task1 = threading.Thread(target = actuateRudder)
task1.start()
fig = plt.figure("press z for zigzag, x for turning, c for free control, a and d controls rudder and press q to quit, space to pause")
# 设置一个网格(grid)，行数为2，列数为2，宽度比例为3:1
gs = gridspec.GridSpec(nrows=2, ncols=2, width_ratios=[1, 1])
ax0 = fig.add_subplot(gs[0, 0])
ax1 = fig.add_subplot(gs[0, 1])
ax2 = fig.add_subplot(gs[1, 0])
ax3 = fig.add_subplot(gs[1, 1])
plt.ion()
while(not stop):
    if(pause):
        ax0.scatter(X, Y,c = 'red')
        ax0.axis('equal')
        ax0.set_title("time:%.1f rudder aim:%.2f deg / rudder now:%.2f deg"%(Time,o_aim,RtoD(data_t[4])))
        ax1.plot(TIME,O_AIM,TIME,THETA)
        ax1.set_title("rudder,heading(deg)")
        ax2.plot(TIME,U,TIME,V)
        ax2.set_title("u,v(normalized)")
        ax3.plot(TIME,R)
        ax3.set_title("angular speed(deg/s)")

        plt.pause(dt)
        ax0.clear()
        ax1.clear()
        ax2.clear()
        ax3.clear()
    else:
        Time = dt*t
        if(zig):
            zigzag(data_t[4],data_t[8])
        if(turn):
            turn10()
        #print(o_aim)
        data_t = realtime_solve(dt,data_t[0],data_t[1],data_t[2],data_t[3],data_t[4],DtoR(o_aim),data_t[6],data_t[7],data_t[8])
        TIME.append(Time)
        X.append(data_t[6]*15)
        Y.append(data_t[7]*15)
        U.append(data_t[0])
        V.append(data_t[1])
        R.append(RtoD(data_t[2]))
        O.append(RtoD(data_t[4]))
        O_AIM.append(o_aim)
        THETA.append(RtoD(data_t[8]))

        ax0.scatter(X, Y,c = 'red')
        ax0.axis('equal')
        ax0.set_title("time:%.1f rudder aim:%.2f deg / rudder now:%.2f deg"%(Time,o_aim,RtoD(data_t[4])))
        ax1.plot(TIME,O_AIM,TIME,THETA)
        ax1.set_title("rudder,heading(deg)")
        ax2.plot(TIME,U,TIME,V)
        ax2.set_title("u,v(normalized)")
        ax3.plot(TIME,R)
        ax3.set_title("angular speed(deg/s)")

        plt.pause(dt)
        ax0.clear()
        ax1.clear()
        ax2.clear()
        ax3.clear()
        t+=1
task1.join()
    
