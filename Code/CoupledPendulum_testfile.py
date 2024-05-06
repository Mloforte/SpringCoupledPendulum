import math
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
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

    def step(self, t,surface):
        # Solve the ODE
        y0 = [self.theta1_0, self.dot_theta1_0,
              self.theta2_0, self.dot_theta2_0]
        sol = odeint(self.G_adim, y0, t)

        # Extract dd_theta1 and dd_theta2
        self.dot_theta1 += sol[:, 1]
        self.dot_theta2 += sol[:, 3]

        self.x1 = self.pivot1[0] + self.l * np.sin(sol[:, 0])
        self.y1 = self.pivot1[1] + self.l * np.cos(sol[:, 0])
        self.x2 = self.pivot2[0] + self.l * np.sin(sol[:, 2])
        self.y2 = self.pivot2[1] + self.l * np.cos(sol[:, 2])

        pygame.draw.line(surface, white, self.pivot1, (self.x1, self.y1))
        pygame.draw.line(surface, white, self.pivot2, (self.x2, self.y2))
        pygame.draw.circle(surface, self.color, (int(
            self.position1[0]), int(self.position1[1])), 15)
        pygame.draw.circle(surface, self.color, (int(
            self.position2[0]), int(self.position2[1])), 15)
        
    def draw(self, surface):
        pygame.draw.line(surface, white, self.pivot1, self.position1)
        pygame.draw.line(surface, white, self.pivot2, self.position2)
        pygame.draw.circle(surface, self.color, (int(
            self.position1[0]), int(self.position1[1])), 15)
        pygame.draw.circle(surface, self.color, (int(
            self.position2[0]), int(self.position2[1])), 15)


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

    while not stop:
        t += 0.1
        clock.tick(fps)
        surface.fill(black)
        for event in pygame.event.get():
            stop = event.type == pygame.QUIT

        pendulum.step()
        pendulum.draw(surface)

        pygame.display.flip()
    pygame.quit()


run()
