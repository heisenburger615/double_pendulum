from vpython import *

canvas(title="Double Pendulum", width=450, height=450, background=color.white, center=vector(0, 0, 0))

# initialize starting conditions
m_1 = 1
m_2 = 1
l_1 = 2
l_2 = 2
g = 9.81
theta_1 = .5 * pi
theta_2 = .75 * pi
d_theta1 = 0
d_theta2 = 0
t = 0
dt = .001



def secDerivs(dtheta1, dtheta2):
    # often used expressions
    sin1 = sin(theta_1)
    sin2 = sin(theta_2)
    cos1 = cos(theta_1)
    cos2 = cos(theta_2)
    sin12 = sin(theta_1 - theta_2)
    cos12 = cos(theta_1 - theta_2)
    m12 = m_1 + m_2

    # differential equations of motion
    A = m12 * l_1
    B = m_2 * l_2 * cos12
    C = -m_2 * l_2 * dtheta2 ** 2 * sin12 - m12 * g * sin1
    D = l_1 / l_2 * cos12
    E = (l_1 * dtheta1 ** 2 * sin12 - g * sin2) / l_2

    dd_theta1 = (C - B * E) / (A - B * D)
    dd_theta2 = E - D * dd_theta1

    energy = -m_1 * g * l_1 * cos1 - m_2 * g * (l_1 * cos1 + l_2 * cos2)
    energy += .5 * m_1 * dtheta1 ** 2 * l_1 ** 2 + .5 * m_2 * (
                dtheta1 ** 2 * l_1 ** 2 + dtheta2 ** 2 * l_2 ** 2 + 2 * l_1 * l_2
                * dtheta1 * dtheta2 * cos12)

    return [dd_theta1, dd_theta2, energy]


# create object
ball1 = sphere(pos=vector(l_1 * sin(theta_1), -l_1 * cos(theta_1), 0),radius=l_1 / 10.0, color=color.blue)
rod1 = cylinder(pos=vector(0, 0, 0), axis=ball1.pos, radius=0.01, color=color.black)
ball2 = sphere(pos=ball1.pos + vector(l_2 * sin(theta_2), -l_2 * cos(theta_2), 0), radius=l_2 / 10.0, color=color.red,
               maketrail=True)
rod2 = cylinder(pos=ball1.pos, axis=ball2.pos - ball1.pos, radius=0.01, color=color.black)

# graph trajectory of second bob
graph(xtitle='x-coordinate', ytitle='y-coordinate')
pos_curve = gcurve(color=color.red)


# motion graphics
while True:
    rate(500)
    dd_theta = secDerivs(d_theta1, d_theta2)  # list of second derivatives

    # Euler-Cromer method
    d_theta1 = d_theta1 + dd_theta[0] * dt
    theta_1 = theta_1 + d_theta1 * dt
    d_theta2 = d_theta2 + dd_theta[1] * dt
    theta_2 = theta_2 + d_theta2 * dt

    t += dt

    ball1.pos = vector(l_1 * sin(theta_1), -l_1 * cos(theta_1), 0)
    rod1.axis = ball1.pos
    ball2.pos = ball1.pos + vector(l_2 * sin(theta_2), -l_2 * cos(theta_2), 0)
    rod2.pos = ball1.pos
    rod2.axis = ball2.pos - ball1.pos

    pos_curve.plot(ball2.pos.x, ball2.pos.y)