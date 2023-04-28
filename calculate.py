import matplotlib.pyplot as plt
import math

# drag constants
p = 1.293
Cd = .4
A = .05
m = .5

# magnus constants
Cl = 0.220 # lift constant

# launch constants
Vxi = 30
Vyi = 10

Vx = Vxi
Vy = Vyi

x = [0]
y = [0]

x1 = [0]
y1 = [0]

def cross(a, b):
    c = [a[1]*b[2] - a[2]*b[1],
         a[2]*b[0] - a[0]*b[2],
         a[0]*b[1] - a[1]*b[0]]

    return c

def get_drag(v):
    # componentes of v
    return -(v/abs(v)) * ((.5*p*(v**2)*Cd*A) / m)

def get_magnus(v, wv, wa):
    """Takes velocity v as tuple of <x,y>, rotational velocity wv, and rotation axis wa as tuple of <x,y>"""
    # vector cross product for direction of magnus force
    # find unit vectors of w and v
    magw = math.sqrt(wa[0] ** 2 + wa[1] ** 2 + wa[2] ** 2)
    magv = math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
    unitv = [v[0] / magv, v[1] / magv, v[2] / magv]
    unitw = [wa[0] / magw, wa[1] / magw, wa[2] / magw]

    magnus_vector = cross(unitv, unitw)
    

    # total v
    return (.5 * Cl * p * A * v**2) / m

dt = .00001
for t in range(0,100000000000000000):
    xa = get_drag(Vx) 
    xy = get_drag(Vy)
    x.append(x[-1] + (Vx * dt + .5 * xa * dt**2))
    y.append(y[-1] + (Vy * dt + .5 * (xy - 9.81) * dt**2))

    x1.append(Vxi * dt * t)
    y1.append(Vyi * dt * t - .5 * 9.8 * (dt*t)**2)

    if y1[-1] < 0:
        break


    Vx += get_drag(Vx) * dt
    Vy += get_drag(Vy) * dt - 9.8 * dt


    if y[-1] < 0:
        break


plt.plot(x,y, 'r')
plt.plot(x1,y1, 'b')
plt.ylim([0, 30])
plt.show()
