import math
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from math import sin, cos, pi, sqrt
import pygame
import sys

black = (0, 0, 0)
white = (255, 255, 255)


class Pendulum:
    def __init__(self, pivot1_x=-10, pivot2_x=10, pivot1_y=0, pivot2_y=0, length=10, mass1=1, mass2=1.5, theta1_0=math.pi/2, theta2_0=math.pi/2, dot_theta1_0=0, dot_theta2_0=0, g=0.981, color='red', k=5.0):
        self.pivot1 = (pivot1_x, pivot1_y)
        self.pivot2 = (pivot2_x, pivot2_y)
        self.l = length
        self.m1 = mass1
        self.m2 = mass2
        self.theta1_0 = theta1_0
        self.theta2_0 = theta2_0
        self.dot_theta1_0 = dot_theta1_0
        self.dot_theta2_0 = dot_theta2_0
        self.g = g
        self.color = color
        self.k = k
        self.x1 = 190
        self.x2 = 220
        self.y1 = 200
        self.y2 = 210

        self.dot_theta1 = 0  # angular velocity
        self.dot_theta2 = 0
        self.trajectory = []
        # Calculation of scaled parameters
        self.A = 3/2
        self.B = (3*self.l*self.k) / (self.g*self.m1)
        self.C = (3*self.l*self.k) / (self.g*self.m2)


    def step(self, t):

        # Solve the ODE
        y0 = [self.theta1_0, self.dot_theta1_0,
              self.theta2_0, self.dot_theta1_0]
        sol = odeint(G_adim, y0, t)

        # Extract dd_theta1 and dd_theta2
        self.dot_theta1 += sol[:, 0]
        sold_theta1 = sol[:, 1]
        self.dot_theta2 += sol[:, 2]
        sold_theta2 = sol[:, 3]

        self.x1 = self.pivot1[0] + self.l * math.sin(self.angle1)
        self.y1 = self.pivot1[1] + self.l * math.cos(self.angle1)
        self.x2 = self.pivot2[0] + self.l * math.sin(self.angle2)
        self.y2 = self.pivot2[1] + self.l * math.cos(self.angle2)

    def draw(self, surface):
        pygame.draw.line(surface, white, self.pivot1, (self.x1, self.y1))
        pygame.draw.line(surface, white, self.pivot2, (self.x2, self.y2))
        pygame.draw.circle(surface, self.color, (self.x1, self.y1), 15)
        pygame.draw.circle(surface, self.color, (self.x2, self.y2), 15)

def G_adim(self,y):
        theta1, dtheta1, theta2, dtheta2 = y
        d_theta1 = dtheta1
        dd_theta1 = self.A*sin(theta1) - self.B*sin(theta1 - theta2) * \
            (1 - (self.d / sqrt(self.d**2 + 2*self.l**2*(1 - cos(theta1 + theta2)))))

        d_theta2 = dtheta2

        dd_theta2 = self.A*sin(theta2) - self.C*sin(theta2 - theta1) * \
            (1 - (self.d / sqrt(self.d**2 + 2*self.l**2*(1 - cos(theta1 + theta2)))))

        return np.array([d_theta1, dd_theta1, d_theta2, dd_theta2])



# # fsolve(g_adim)

# # Solve the ODE
# y0 = [theta1_0, dot_theta1_0, theta2_0, dot_theta1_0]
# sol = odeint(G_adim, y0, t)

# # Initialize plot
# plt.ion()
# fig, ax = plt.subplots()
# line1, = ax.plot([], [], label='theta1')
# line2, = ax.plot([], [], label='theta2')
# ax.set_xlabel('Temps')
# ax.set_ylabel('Theta')
# ax.set_title('Oscillations du système en fonction du temps')
# ax.legend()
# ax.grid(True)

# # Update plot in a loop
# t = 0
# while True:
#     # Increment time
#     t += 0.1

#     # Solve the ODE
#     y0 = [theta1_0, dot_theta1_0, theta2_0, dot_theta2_0]
#     sol = odeint(G_adim, y0, [0, t])

#     # Extract theta1 and theta2
#     theta1 = sol[-1, 0]
#     theta2 = sol[-1, 2]

#     # Update plot data
#     line1.set_xdata(np.append(line1.get_xdata(), t))
#     line1.set_ydata(np.append(line1.get_ydata(), theta1))
#     line2.set_xdata(np.append(line2.get_xdata(), t))
#     line2.set_ydata(np.append(line2.get_ydata(), theta2))

#     # Redraw the plot
#     ax.relim()
#     ax.autoscale_view()

#     # Pause to give time to redraw
#     plt.pause(0.01)

#     # Guess initial values for the stationary points
#     # You can provide initial guesses based on your system
#     y_guess = [0, 0, 0, 0]

#     # Solve for the stationary points
#     stationary_solution = fsolve(G_adim, y_guess, t)

#     # Extract the stationary points
#     theta1_stationary, dtheta1_stationary, theta2_stationary, dtheta2_stationary = stationary_solution


# pendulum loop
def init_surface(size, caption):
    pygame.init()
    pygame.display.set_caption(caption)
    surface = pygame.display.set_mode(size)
    clock = pygame.time.Clock()
    return surface, clock


def run():
    width, height = 800, 800
    fps = 30
    surface, clock = init_surface((width, height), 'Double Pendulum')

    pendulum = Pendulum(width//2, width//2, height//2, height//2)
    stop = False
    # Time points to solve the ODE for
    t = np.linspace(0, 10, 1000)

    while not stop:
        t += 0.1
        clock.tick(fps)
        surface.fill(black)
        for event in pygame.event.get():
            stop = event.type == pygame.QUIT

        pendulum.step(t)
        pendulum.draw(surface)
        pygame.display.flip()
    pygame.quit()


run()
