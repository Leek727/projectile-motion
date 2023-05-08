import matplotlib.pyplot as plt
import math
import numpy as np
# drag constants
p = 1.293
Cd = .4
R = .1 # radius
A = math.pi * R ** 2
m = .270

# magnus constants
Cl = 0.220  # lift constant

# launch constants
Vxi = 5
Vyi = 5
Vzi = 1
w = 1

# used vars
Vx = Vxi
Vy = Vyi
Vz = Vzi

x0 = [0]
y0 = [0]
z0 = [0]

x1 = [0]
y1 = [0]
z1 = [0]

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


# current position for drag
x = 0
y = 0
z = 0

drag_done = False
vac_done = False

dt = .001
t = 0
print("Starting... ")
while True:
    # magnus and drag effect numerical integration TODO use better integration method
    # x,y,z
    am = get_magnus([Vz, Vy, Vx], 4000, [1, 1, 0])  # [1,0,0])

    ax = get_drag(Vx) + am[2]
    ay = get_drag(Vy) - 9.81 + am[1]
    az = get_drag(Vz) + am[0]

    Vxold = Vx
    Vyold = Vy
    Vzold = Vz

    Vx += ax * dt
    Vy += ay * dt
    Vz += az * dt

    x += dt * (Vxold + Vx) / 2 # i like trapez
    y += dt * (Vyold + Vy) / 2
    z += dt * (Vzold + Vz) / 2

    if y0[-1] < 0 and not drag_done:
        print(f"Ball with drag lands at: {round(t * dt,3)} s")
        drag_done = True
    if not drag_done:
        x0.append(x)
        y0.append(y)
        z0.append(z)

    if y1[-1] < 0 and not vac_done:
        print(f"Ball without drag lands at: {round(t * dt,3)} s")
        vac_done = True

    if not vac_done:
        # vacuum
        x1.append(Vxi * dt * t)
        y1.append(Vyi * dt * t - .5 * 9.8 * (dt*t)**2)
        z1.append(Vzi * dt * t)

    if y1[-1] < 0 and y0[-1] < 0:
        break

    t += 1

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


ax = plt.axes(projection ='3d')
ax.plot3D(x0, z0, y0, 'red')
ax.plot3D(x1, z1, y1, 'blue')

set_axes_equal(ax)
plt.show()
