import matplotlib.pyplot as plt
import math

# drag constants
p = 1.293
Cd = .4
A = 0.045
m = .270

# magnus constants
Cl = 0.220  # lift constant

# launch constants
angle = 60
Vxi = 17.3736 * math.cos((angle/180) * math.pi)
Vyi = 17.3736 * math.sin((angle/180) * math.pi)
print(Vxi)
w = 1

# used vars
Vx = Vxi
Vy = Vyi

x0 = [0]
y0 = [0]

x1 = [0]
y1 = [0]


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

    # S = wv/magv # spin parameter.... http://spiff.rit.edu/richmond/baseball/traj/traj.html apparently affects Cl somehow experimentally???????/

    # total v
    return am


# current position for drag
x = 0
y = 0

drag_done = False
vac_done = False

dt = .001
t = 0
while True:
    # magnus and drag effect numerical integration TODO use better integration method
    # x,y,z - TODO make entire sim 3d
    am = get_magnus([0, Vy, Vx], 1, [1, 0, 0])  # [1,0,0])

    ax = get_drag(Vx) + am[2]
    ay = get_drag(Vy) - 9.81 + am[1]

    Vxold = Vx
    Vyold = Vy

    Vx += ax * dt
    Vy += ay * dt

    x += dt * (Vxold + Vx) / 2  # trap
    y += dt * (Vyold + Vy) / 2

    if y0[-1] < 0 and not drag_done:
        print(f"Ball with drag lands at: {round(t * dt,3)} s")
        drag_done = True
    if not drag_done:
        x0.append(x)
        y0.append(y)

    if y1[-1] < 0 and not vac_done:
        print(f"Ball without drag lands at: {round(t * dt,3)} s")
        vac_done = True

    if not vac_done:
        # vacuum
        x1.append(Vxi * dt * t)
        y1.append(Vyi * dt * t - .5 * 9.8 * (dt*t)**2)

    if y1[-1] < 0 and y0[-1] < 0:
        break

    t += 1

plt.plot(x0, y0, 'r')
plt.plot(x1, y1, 'b')
plt.axis('scaled')
#plt.ylim([0, 30])
plt.show()
