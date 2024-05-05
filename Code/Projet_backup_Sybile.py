import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from math import sin, cos, pi, sqrt
import pygame
import sys

# Initialisation de Pygame
pygame.init()

# Création de la fenêtre
screen = pygame.display.set_mode((400, 300))
pygame.display.set_caption('Simulation du système pendulaire')

Out = False  # if True,out of while loop, and close pygame
# System constants
m1 = m2 = 1.5  # masse du pendule
l = 15  # longueur de la suspension
d = 30  # distance entre les pendules où le ressort est attaché
k = 5.0  # coefficient de raideur du ressort
g = 98  # accélération due à la gravité


# Initial conditions
theta1_0 = pi/2  # déviation initiale du premier pendule
theta2_0 = pi/2  # déviation initiale du deuxième pendule
dot_theta1_0 = 0.0  # vitesse angulaire initiale du premier pendule
dot_theta2_0 = 0.0  # vitesse angulaire initiale du deuxième pendule

# Calculation of scaled parameters
A = 3/2
B = (3*l*k) / (g*m1)
C = (3*l*k) / (g*m2)

# equation =

# solution = sol(equation, x)


def G_adim(y, t):
    theta1, dtheta1, theta2, dtheta2 = y
    d_theta1 = dtheta1
    dd_theta1 = A*sin(theta1) - B*sin(theta1 - theta2) * \
        (1 - (d / sqrt(d**2 + 2*l**2*(1 - cos(theta1 + theta2)))))
    d_theta2 = dtheta2
    dd_theta2 = A*sin(theta2) - C*sin(theta2 - theta1) * \
        (1 - (d / sqrt(d**2 + 2*l**2*(1 - cos(theta1 + theta2)))))
    return np.array([d_theta1, dd_theta1, d_theta2, dd_theta2])


# Time points to solve the ODE for
t = np.linspace(0, 10, 1000)

# fsolve(g_adim)

# Solve the ODE
y0 = [theta1_0, dot_theta1_0, theta2_0, dot_theta1_0]
sol = odeint(G_adim, y0, t)

# Extract dd_theta1 and dd_theta2
soltheta1 = sol[:, 0]
sold_theta1 = sol[:, 1]
soltheta2 = sol[:, 2]
sold_theta2 = sol[:, 3]

# Plot dd_theta1 and dd_theta2
# plt.plot(t,soltheta1, label='theta1')
# plt.plot(t, soltheta2, label='theta2')
# plt.xlabel('Temps')
# plt.ylabel('theta')
# plt.title('theta1 et theta2 en fonction du temps')
# plt.legend()
# plt.grid(True)
# plt.show()

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

# Update plot in a loop
t = 0
while not Out:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            Out = True  # Si la fenêtre est fermée, sortir de la boucle
    # Increment time
    t += 0.1

    # Solve the ODE
    y0 = [theta1_0, dot_theta1_0, theta2_0, dot_theta2_0]
    sol = odeint(G_adim, y0, [0, t])

    # Extract theta1 and theta2
    theta1 = sol[-1, 0]
    theta2 = sol[-1, 2]

    # Update plot data
    line1.set_xdata(np.append(line1.get_xdata(), t))
    line1.set_ydata(np.append(line1.get_ydata(), theta1))
    line2.set_xdata(np.append(line2.get_xdata(), t))
    line2.set_ydata(np.append(line2.get_ydata(), theta2))

    # Redraw the plot
    ax.relim()
    ax.autoscale_view()

    # Pause to give time to redraw
    plt.pause(0.01)

    # Guess initial values for the stationary points
    # You can provide initial guesses based on your system
    y_guess = [0, 0, 0, 0]

    # Solve for the stationary points
    stationary_solution = fsolve(G_adim, y_guess, t)

    # Extract the stationary points
    theta1_stationary, dtheta1_stationary, theta2_stationary, dtheta2_stationary = stationary_solution

pygame.quit()
