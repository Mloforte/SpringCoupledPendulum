# Wesley Fernandes
# python simple pendulum with pygame

import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import odeint
import pygame
import math
from math import sin, cos, pi, sqrt
import matplotlib.pyplot as plt

# VARIABLES
width, height = 800, 400   # set the width and height of the window
# (you can increase or decrease if you want to, just remind to keep even numbers)
Out = False                # if True,out of while loop, and close pygame
# when true it allow us to find the acceleration and damping for the pendulum
acceleration = False
length = 15                # the length between the ball and the support
angle1 = pi/4                  # the angle that you begin when click in window
vel1 = 0                    # velocity that angle is increased and damped
Aacc1 = 0                   # acceleration
m1 = 1.5
angle2 = pi/2                  # the angle that you begin when click in window
vel2 = 0                    # velocity that angle is increased and damped
Aacc2 = 0                   # acceleration
m2 = 1.5
g = 9.81
k = 5.0
d = 30
A = 3/2
B = (3*length*k) / (g*m1)
C = (3*length*k) / (g*m2)

# COLORS
white = (255, 255, 255)
black = (0, 0, 0)
gray = (150, 150, 150)
Dark_red = (150, 0, 0)

# BEFORE START
pygame.init()
background = pygame.display.set_mode((width, height))
clock = pygame.time.Clock()


class ball(object):

    def __init__(self, XY1, XY2, radius):  # Set ball coordenates and radius
        self.x1 = XY1[0]
        self.y1 = XY1[1]
        self.x2 = XY2[0]
        self.y2 = XY2[1]
        self.radius = radius

    def draw(self, bg):  # Draw circle and line based on XY coordinates
        pygame.draw.lines(bg, black, False, [
                          (-d/2 + 400, 50), (self.x1, self.y1)], 2)
        pygame.draw.circle(bg, black, (self.x1, self.y1), self.radius)
        # pygame.draw.circle(bg, Dark_red, (self.x1, self.y1), self.radius - 2)
        pygame.draw.lines(bg, black, False, [
                          (d/2 + 400, 50), (self.x2, self.y2)], 2)
        pygame.draw.circle(bg, black, (self.x2, self.y2), self.radius)
        # pygame.draw.circle(bg, Dark_red, (self.x2, self.y2), self.radius - 2)
        pygame.draw.lines(bg, gray, False, [
                          (self.x1, self.y1), (self.x2, self.y2)], 2)


def calclPosBall2(XY1, d):
    position = (XY1[0] + d, XY1[1])
    return position


def grid():  # Draw a grid behind the pendulum
    for x in range(50, width, 50):
        pygame.draw.lines(background, gray, False, [(x, 0), (x, height)])
        for y in range(50, height, 50):
            pygame.draw.lines(background, gray, False, [(0, y), (width, y)])
    # pygame.draw.circle(background, black, (int(width/2), 50), 5)


def angle_Length():  # Send back the length and angle at the first click on screen
    # length = math.sqrt(math.pow(pendulum.x - width/2, 2) +math.pow(pendulum.y - 50, 2))

    angle = math.asin((pendulum.x - width/2) / length)
    return (angle, length)


def get_path(angle1, angle2, length):  # with angle and length calculate x and y position
    pendulum.x1 = round(400 - (d/2) + (length * math.sin(angle1)))
    pendulum.y1 = round(50 + (length* math.cos(angle1)))
    pendulum.x2 = round( 400+ (d/2) + (length * math.sin(angle2)))
    pendulum.y2 = round(50 + (length * math.cos(angle2)))


def redraw():  # Clean up the screen and start a new grid and new frame of pendulum with new coordinates
    background.fill(white)
    grid()
    pendulum.draw(background)
    pygame.display.update()


def G_adim(y, t):

    angle1, dangle1, angle2, dangle2 = y
    d_angle1 = dangle1
    dd_angle1 = A * sin(angle1) - B * sin(angle1 - angle2) * \
        (1 - (d / sqrt(d**2 + 2 * length**2 * (1 - cos(angle1 + angle2)))))

    d_angle2 = dangle2
    dd_angle2 = A * sin(angle2) - B * sin(angle2 - angle1) * \
        (1 - (d / sqrt(d**2 + 2 * length**2 * (1 - cos(angle1 + angle2)))))

    return np.array([d_angle1, dd_angle1, d_angle2, dd_angle2])


# I start the class with some random coordinates
pendulum = ball((-d/2 + 400, length), (d/2 + 400, length), 5)

# Time points to solve the ODE for
t = np.linspace(0, 10, 1000)
t = 0
# Initialize plot
plt.ion()
fig, ax = plt.subplots()
line1, = ax.plot([], [], label='theta1')
line2, = ax.plot([], [], label='theta2')
ax.set_xlabel('Temps')
ax.set_ylabel('Theta')
ax.set_title('Oscillations du système en fonction du temps')
ax.legend()
ax.grid(True)

while not Out:
    clock.tick(120)  # Set how many frames are draw per second
    # If changed, maybe, could be a good idea change some values at acceleration

    # Increment time
    t += 0.1

    for event in pygame.event.get():                     #
        if event.type == pygame.QUIT:                    #
            Out = True                                   #
        if event.type == pygame.MOUSEBUTTONDOWN:  # Read if you want go out
            # pendulum = ball(pygame.mouse.get_pos(), 15)  # or
            pendulum = ball((-50, -50), calclPosBall2((-50, -50), d), 5)
            # angle, length = angle_Length()  # click the mouse button
            # angles and length defined above
            acceleration = True

    if acceleration:   # Increase acceleration and damping in the pendulum moviment

        # Aacc1 = -0.005 * math.sin(angle1)
        # vel1 += Aacc1
        # vel1 *= 0.99  # damping factor
        # angle1 += vel1

        # Aacc2 = -0.005 * math.sin(angle2)
        # vel2 += Aacc2
        # vel2 *= 0.99  # damping factor
        # angle2 += vel2

        # Solve the ODE
        y0 = [angle1, vel1, angle2, vel2]
        sol = odeint(G_adim, y0, [0, t])
        # Extract angle1, vel1, angle2, vel2 from the solution
        angle1 = sol[-1, 0]
        vel1 = sol[-1, 1]
        angle2 = sol[-1, 2]
        vel2 = sol[-1, 3]
        get_path(angle1, angle2, length)
        line1.set_xdata(np.append(line1.get_xdata(), t))
        line1.set_ydata(np.append(line1.get_ydata(), angle1))
        line2.set_xdata(np.append(line2.get_xdata(), t))
        line2.set_ydata(np.append(line2.get_ydata(), angle2))

        # Redraw the plot
        ax.relim()
        ax.autoscale_view()

        # Pause to give time to redraw
        plt.pause(0.01)

    redraw()

pygame.quit()
