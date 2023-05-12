import matplotlib.pyplot as plt
import math
import numpy as np
import time
import random

# drag constants
p = 1.293
Cd = .4
R = .1 # radius
A = math.pi * R ** 2
m = .270

# magnus constants
Cl = 0.220  # lift constant


def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])


def cross(a, b):
    """Takes cross product of two vectors and returns a unit vector"""
    c = [a[1]*b[2] - a[2]*b[1],
         a[2]*b[0] - a[0]*b[2],
         a[0]*b[1] - a[1]*b[0]]

    return c


def get_drag(v):
    # componentes of v
    return -(v/abs(v)) * ((.5*p*(v**2)*Cd*A) / m)


def get_magnus(v, wv, wa):
    if wv == 0:
        return [0,0,0]
    """Takes velocity v as tuple of <x,y,z>, rotational velocity wv, and rotation axis wa as tuple of <x,y,z>"""
    # vector cross product for direction of magnus force
    magw = math.sqrt(wa[0] ** 2 + wa[1] ** 2 + wa[2] ** 2)
    magv = math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
    unitv = [v[0] / magv, v[1] / magv, v[2] / magv]
    unitw = [wa[0] / magw, wa[1] / magw, wa[2] / magw]

    magnus_hat = cross(unitv, unitw)  # magnus carlsen vector wearing a hat

    # find Fm and orient
    vmag = math.sqrt(v[0]**2 + v[1]**2)
    am = (.5 * Cl * p * A * vmag**2) / m
    am = [am * magnus_hat[0], am * magnus_hat[1], am *
          magnus_hat[2]]  # orient towards unit vector

    # EXPERIMENTAL STUF
    S = (R * wv)/magv # spin parameter.... http://spiff.rit.edu/richmond/baseball/traj/traj.html apparently affects Cl somehow experimentally???????/

    # total v
    return [S * x for x in am]



def get_acc(Vx, Vy, Vz, w, unit_direction):
    # x,y,z
    am = get_magnus([Vz, Vy, Vx], w, unit_direction)  # [1,0,0])

    ax = get_drag(Vx) + am[2]
    ay = get_drag(Vy) - 9.81 + am[1]
    az = get_drag(Vz) + am[0]
    return [ax, ay, az]


def rungekutta_step(f, Vx, Vy, Vz, dt, w, unit_direction):
    """Take as input current velocities"""
    start_val = np.array([Vx, Vy, Vz])

    k1 = dt * np.array(f(Vx, Vy, Vz, w, unit_direction))
    
    a,b,c = start_val + k1/2
    k2 = dt * np.array(f(a,b,c, w, unit_direction))
    
    a,b,c = start_val + k2 / 2
    k3 = dt * np.array(f(a,b,c, w, unit_direction))
    
    a,b,c = start_val + k3
    k4 = dt * np.array(f(a,b,c, w, unit_direction))
    final = start_val + (k1 + 2 * k2 + 2 * k3 + k4) / 6
    
    return final


def velocity_step(f, Vx, Vy, Vz, dt, w, unit_direction):
    """Take as input current velocities"""
    start_val = np.array([Vx, Vy, Vz])

    k1 = dt * np.array(f(Vx, Vy, Vz, w, unit_direction))
    
    a,b,c = start_val + k1/2
    k2 = dt * np.array(f(a,b,c, w, unit_direction))
    
    a,b,c = start_val + k2 / 2
    k3 = dt * np.array(f(a,b,c, w, unit_direction))
    
    a,b,c = start_val + k3
    k4 = dt * np.array(f(a,b,c, w, unit_direction))
    final = start_val + (k1 + 2 * k2 + 2 * k3 + k4) / 6
    
    return final


def calc(Vx, Vy, Vz, w, unit_direction, runge=False):
    Vx, Vy, Vz, w = Vx + .00000001, Vy + .00000001, Vz + .00000001, w + .00000001
    x0 = [0]
    y0 = [0]
    z0 = [0]

    # current position for drag
    x = 0
    y = 0
    z = 0

    drag_done = False

    dt = .1
    t = 0
    while True:

        Vxold = Vx
        Vyold = Vy
        Vzold = Vz

        if runge:
            dt = .2
            Vx, Vy, Vz = rungekutta_step(get_acc, Vx, Vy, Vz, dt, w, unit_direction)
        
        else:
            ax,ay,az = get_acc(Vx, Vy, Vz, w, unit_direction)

            Vx += ax * dt
            Vy += ay * dt
            Vz += az * dt
    


        #print(np.array([tx,ty,tz]) - np.array([Vx, Vy, Vz]))

        x += dt * (Vxold + Vx) / 2 # i like trapez
        y += dt * (Vyold + Vy) / 2
        z += dt * (Vzold + Vz) / 2

        
        # check target
        # z0 depth, y0 height, x0 horizontal
        """  if z0[-1] > 5:
            if y0[-1] > 2 and y0[-1] < 2.5 and x0[-1] > -.5 and x0[-1] < .5:
                return [x0, y0, z0]

            break"""
        


        if y0[-1] < 0 and not drag_done:
            #print(f"Ball with drag lands at: {round(t * dt,3)} s")
            return [x0, y0, z0]
            break

        if not drag_done:
            x0.append(x)
            y0.append(y)
            z0.append(z)

        t += 1

    return [0,0,0]

ax = plt.axes(projection ='3d')

"""
perms = [
        [0,0,1],
        [0,1,0],
        [0,1,1],
        [1,0,0],
        [1,0,1],
        [1,1,1],    
]
# omni source
for x in range(-20, 20, 2):
    y = 5
    for z in range(-20, 20, 2):
        for j in range(-5,5):
            for perm in perms:
                x0, y0, z0 = calc(x, y, z, j * 400, perm)
                if x0 + y0 + z0 != 0:
                    ax.plot3D(x0, z0, y0, 'blue')

    print(x)
"""



"""for i in range(-2,3):
    for j in range(-2,3):
        for k in range(-2,3):
            mag = [i+.00000001,j+.00000001,k+.00000001]
            x0, y0, z0 = calc(10,10,0, 100, mag, False)
            if x0 + y0 + z0 != 0:
                ax.plot3D(x0, z0, y0, 'blue')

for i in range(-2,3):
    for j in range(-2,3):
        for k in range(-2,3):
            mag = [i+.00000001,j+.00000001,k+.00000001]
            x0, y0, z0 = calc(10,10,0, 100, mag, True)
            if x0 + y0 + z0 != 0:
                ax.plot3D(x0, z0, y0, 'red')

    print(i)
"""

ground_truth = np.array([14.361254594177511, -3.675186128917694e-05, 2.862851170238706])

s = time.time()
x0, y0, z0 = calc(10,10,0, 1000, [0,0,1], False)
if x0 + y0 + z0 != 0:
    ax.plot3D(x0, z0, y0, 'blue')

e = time.time() - s
print(f"Error : {sum([abs(x) for x in ground_truth - np.array([x0[-1],y0[-1],z0[-1]])])}")
print(f"Trapezoidal method, t = {e}")


s = time.time()
x0, y0, z0 = calc(10,10,0, 1000, [0,0,1], True)
if x0 + y0 + z0 != 0:
    ax.plot3D(x0, z0, y0, 'red')

e = time.time() - s
print(f"Error : {sum([abs(x) for x in ground_truth - np.array([x0[-1],y0[-1],z0[-1]])])}")
print(f"Runge Kutta method, t = {e}")


set_axes_equal(ax)
plt.show()
