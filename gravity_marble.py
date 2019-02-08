import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation, rc
from IPython.display import HTML


matplotlib.rcParams['animation.writer'] = 'avconv'

plt.rcParams['figure.figsize'] = (6,3)
pipex = np.linspace(-1,1,100)
pipey = -(1-(pipex**2))**0.5
plt.plot(pipex,pipey)

def rk4(f,t,h,g,**kwargs):
    k1 = h*g(t,f,**kwargs)
    k2 = h*g(t+0.5*h, f+0.5*k1,**kwargs)
    k3 = h*g(t+0.5*h, f+0.5*k2,**kwargs)
    k4 = h*g(t+h, f+k3,**kwargs)
    
    return f+ k1/6. + k2/3. + k3/3. + k4/6.

def CircularHarmOsc(t,f):
    x = f[0]
    v = f[1]
    gravity = 9.81
    RHS = np.array([v,-gravity*x])
    return RHS

# initial conditions:
x0 = np.array([-1.,0.])
t0 = 0. # start at t= 0

t_stop = 15
x = [x0[0],]
v = [x0[1],]
t = [t0,]
xold = x0
h=0.1
while t[-1] < t_stop:
    xold = rk4(xold,t[-1],h,CircularHarmOsc)
    x.append(xold[0])
    v.append(xold[1])
    t.append(t[-1]+h)
x = np.array(x)
t = np.array(t)
v = np.array(v)

plt.figure()
plt.plot(t,x, label="Position")
plt.plot(t,v, label="Velocity")
plt.legend(loc = 'upper right')
plt.xlabel("t")
plt.ylabel("x or v")

fig, ax = plt.subplots()

ax.plot(pipex,pipey)
board = ax.scatter(-1,0)

def animate(i):
    board.set_offsets([x[i],y[i]])


y = -(1-(x**2))**0.5

anim = animation.FuncAnimation(
    fig, animate, interval=100, frames=len(t)-1)

HTML(anim.to_html5_video())

def DampCircularHarmOsc(t,f,beta = 0.1,gamma = 0.4):
    x = f[0]
    v = f[1]
    gravity = 9.81
    
    RHS = np.array([v,-gravity*x -(beta + gamma)*v])
    return RHS

# initial conditions:
x0 = np.array([-1.,0.])
t0 = 0. # start at t= 0
beta = 0.1
gamma = 0.4
t_stop = 15
x = [x0[0],]
v = [x0[1],]
t = [t0,]
xold = x0
h=0.1
while t[-1] < t_stop:
    xold = rk4(xold,t[-1],h,DampCircularHarmOsc,beta = beta, gamma = gamma)
    x.append(xold[0])
    v.append(xold[1])
    t.append(t[-1]+h)
x = np.array(x)
v = np.array(v)
t = np.array(t)

plt.figure()
plt.plot(t,x, label="Position")
plt.plot(t,v, label="Velocity")
plt.legend(loc = 'upper right')
plt.xlabel("t")
plt.ylabel("x or v")


fig, ax = plt.subplots()
ax.plot(pipex,pipey)
board = ax.scatter(-1,0)
y = -(1-(x**2))**0.5
anim = animation.FuncAnimation(
    fig, animate, interval=100, frames=len(t)-1)

HTML(anim.to_html5_video())

def CircularPipe(d):
    CircularX = np.linspace(-1,1,100)
    CircularY = -(((d**2 + 1)/(2*d))**2-(CircularX**2))**0.5 + ((1-(d**2))/(2*d))
    plt.plot(CircularX,CircularY)
    plt.xlim(-1.1,1.1)
    plt.ylim(-1.1,0.1)

CircularPipe(0.6)

def GeneralCircularHarmOsc(t,f,beta = 0.1,gamma = 0.4,d = 0.6):
    x = f[0]
    v = f[1]
    gravity = 9.81
    
    RHS = np.array([v,-gravity*((2*d*x)/((d**2)+1)) -(beta + gamma)*v])
    return RHS

x0 = np.array([-1.,0.])
t0 = 0. # start at t= 0
beta = 0.1
gamma = 0.4
d = 0.6
t_stop = 15
x = [x0[0],]
v = [x0[1],]
t = [t0,]
xold = x0
h=0.1
while t[-1] < t_stop:
    xold = rk4(xold,t[-1],h,GeneralCircularHarmOsc,beta = beta, gamma = gamma,d = d)
    x.append(xold[0])
    v.append(xold[1])
    t.append(t[-1]+h)
x = np.array(x)
v = np.array(v)
t = np.array(t)

plt.figure()
plt.plot(t,x, label="Position")
plt.plot(t,v, label="Velocity")
plt.legend(loc = 'upper right')
plt.xlabel("t")
plt.ylabel("x or v")

CircularX = np.linspace(-1,1,100)
CircularY = -(((0.6**2 + 1)/(2*0.6))**2-(CircularX**2))**0.5 + ((1-(0.6**2))/(2*0.6))

fig, ax = plt.subplots()
ax.plot(CircularX,CircularY)
plt.xlim(-1.1,1.1)
plt.ylim(-1.1,0.1)
board = ax.scatter(-1,0)
y = -(((0.6**2 + 1)/(2*0.6))**2-(x**2))**0.5 + ((1-(0.6**2))/(2*0.6))
anim = animation.FuncAnimation(
    fig, animate, interval=100, frames=len(t)-1)

HTML(anim.to_html5_video())

EOverMCircular = 9.81*(y+0.6) + 0.5*(v**2)

SinePipeX = np.linspace(-1,1,100)
SinePipeY = 0.5*np.cos(np.pi*SinePipeX + np.pi) - 0.5
plt.plot(SinePipeX,SinePipeY)

def SineHarmOsc(t,f):
    x = f[0]
    v = f[1]
    gravity = 9.81
    dydx = -0.5*np.pi*np.sin(np.pi*x + np.pi)
    RHS = np.array([v,-gravity*np.sin(np.arctan(dydx))])
    return RHS

# initial conditions:
x0 = np.array([-1.,0.0001])
t0 = 0. # start at t= 0

t_stop = 10
x = [x0[0],]
v = [x0[1],]
t = [t0,]
xold = x0
h=0.1
while t[-1] < t_stop:
    xold = rk4(xold,t[-1],h,SineHarmOsc)
    x.append(xold[0])
    v.append(xold[1])
    t.append(t[-1]+h)
x = np.array(x)
v = np.array(v)
t = np.array(t)

plt.figure()
plt.plot(t,x, label="Position")
plt.plot(t,v, label="Velocity")
plt.legend(loc = 'upper right')
plt.xlabel("t")
plt.ylabel("x or v")

fig, ax = plt.subplots()
ax.plot(SinePipeX,SinePipeY)
board = ax.scatter(-1,0)
y = 0.5*np.cos(np.pi*x + np.pi) - 0.5
anim = animation.FuncAnimation(
    fig, animate, interval=100, frames=len(t)-1)

HTML(anim.to_html5_video())

def SineDampHarmOsc(t,f,beta = 0.1,gamma = 0.4):
    x = f[0]
    v = f[1]
    gravity = 9.81
    dydx = -0.5*np.pi*np.sin(np.pi*x + np.pi)
    RHS = np.array([v,-gravity*np.sin(np.arctan(dydx))-(beta + gamma)*v])
    return RHS

# initial conditions:
x0 = np.array([-1.,0.0001])
t0 = 0. # start at t= 0
beta = 0.1
gamma = 0.4
t_stop = 10
x = [x0[0],]
v = [x0[1],]
t = [t0,]
xold = x0
h=0.1
while t[-1] < t_stop:
    xold = rk4(xold,t[-1],h,SineDampHarmOsc, beta = beta, gamma = gamma)
    x.append(xold[0])
    v.append(xold[1])
    t.append(t[-1]+h)
x = np.array(x)
t = np.array(t)
v = np.array(v)

plt.figure()
plt.plot(t,x, label="Position")
plt.plot(t,v, label="Velocity")
#plt.axhline(0)
plt.legend(loc = 'upper right')
plt.xlabel("t")
plt.ylabel("x or v")

fig, ax = plt.subplots()
ax.plot(SinePipeX,SinePipeY)
board = ax.scatter(-1,0)
y = 0.5*np.cos(np.pi*x + np.pi) - 0.5
anim = animation.FuncAnimation(
    fig, animate, interval=100, frames=len(t)-1)

HTML(anim.to_html5_video())

def SinusoidalPipe(d):
    SineX = np.linspace(-1,1,100)
    SineY = (d/2)*np.cos(np.pi*SineX + np.pi) - (d/2)
    plt.plot(SineX,SineY)
    plt.xlim(-1.1,1.1)
    plt.ylim(-1.1,0.1)

SinusoidalPipe(0.6)

def GeneralSineHarmOsc(t,f,beta = 0.1,gamma = 0.4,d = 0.6):
    x = f[0]
    v = f[1]
    gravity = 9.81
    dydx = -(d/2)*np.pi*np.sin(np.pi*x + np.pi)
    RHS = np.array([v,-gravity*np.sin(np.arctan(dydx)) -(beta + gamma)*v])
    return RHS

x0 = np.array([-1.,0.0001])
t0 = 0. # start at t= 0
beta = 0.1
gamma = 0.4
d = 0.6
t_stop = 15
x = [x0[0],]
v = [x0[1],]
t = [t0,]
xold = x0
h=0.1
while t[-1] < t_stop:
    xold = rk4(xold,t[-1],h,GeneralSineHarmOsc,beta = beta, gamma = gamma,d = d)
    x.append(xold[0])
    v.append(xold[1])
    t.append(t[-1]+h)
x = np.array(x)
v = np.array(v)
t = np.array(t)

plt.figure()
plt.plot(t,x, label="Position")
plt.plot(t,v, label="Velocity")
plt.legend(loc = 'upper right')
plt.xlabel("t")
plt.ylabel("x or v")

SineX = np.linspace(-1,1,100)
SineY = (0.6/2)*np.cos(np.pi*SineX + np.pi) - (0.6/2)

fig, ax = plt.subplots()
ax.plot(SineX,SineY)
plt.xlim(-1.1,1.1)
plt.ylim(-1.1,0.1)
board = ax.scatter(-1,0)
y = (0.6/2)*np.cos(np.pi*x + np.pi) - (0.6/2)
anim = animation.FuncAnimation(
    fig, animate, interval=100, frames=len(t)-1)

HTML(anim.to_html5_video())

EOverMSinusoidal = 9.81*(y+0.6) + 0.5*(v**2)
plt.plot(t,EOverMCircular, label="CircularPipe")
plt.plot(t,EOverMSinusoidal, label="SinePipe")
plt.xlabel("time")
plt.ylabel("Total energy/mass")
plt.legend(loc = 'upper right')
