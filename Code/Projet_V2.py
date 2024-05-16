import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from math import sin, cos, pi, sqrt
import pygame
import sys
from scipy.fft import fft, fftfreq


width, height = 800, 400

# # Initialisation of Pygame
pygame.init()

# # Creation of the fenêtre
background = pygame.display.set_mode((800, 800))

# COLORS
white = (255, 255, 255)
black = (0, 0, 0)
gray = (150, 150, 150)
Dark_red = (150, 0, 0)

# Loop parameters
Out = False  # if True,out of while loop, and close pygame

# System constants
m1 = 100  # pendulum mass 1
m2 = 100  # pendulum mass 2
l = 50  # suspension length
d = 100  # distance between pendulums where the spring is attached
k = 5.0  # coefficient of spring stiffness
g = 9.8  # acceleration due to gravity
x = 0.5  # damping factor for lose of energy due to friction
acceleration = False

# Initial conditions
theta1_0 = -pi/2  # deviation of the first pendulum
theta2_0 = -pi/4  # deviation of the second pendulum
dot_theta1_0 = 0.0  # intial angular velocity of the first pendulum
dot_theta2_0 = 0.0  # initial angular velocity of the second pendulum

# Calculation of scaled parameters
A = 3/2
B = (3*l*k) / (g*m1)
C = (3*l*k) / (g*m2)


class ball(object):

    def __init__(self, XY1, XY2, radius):  # Set ball coordenates and radius
        self.x1 = XY1[0]
        self.y1 = XY1[1]
        self.x2 = XY2[0]
        self.y2 = XY2[1]
        self.radius = radius

    def draw(self, bg):  # Draw circle and line based on XY coordinates
        pygame.draw.lines(bg, black, False, [
                          (-d/2 + width/2, height/4), (self.x1, self.y1)], 2)
        pygame.draw.circle(bg, black, (self.x1, self.y1), self.radius)
        # pygame.draw.circle(bg, Dark_red, (self.x1, self.y1), self.radius - 2)
        pygame.draw.lines(bg, black, False, [
                          (d/2 + width/2, height/4), (self.x2, self.y2)], 2)
        pygame.draw.circle(bg, black, (self.x2, self.y2), self.radius)
        # pygame.draw.circle(bg, Dark_red, (self.x2, self.y2), self.radius - 2)
        pygame.draw.lines(bg, gray, False, [
                          (self.x1, self.y1), (self.x2, self.y2)], 2)


def get_path(angle1, angle2, length):  # with angle and length calculate x and y position
    pendulum.x1 = round(width/2 - (d/2) + (length * sin(angle1)))
    pendulum.y1 = round(height/4 + (length * cos(angle1)))
    pendulum.x2 = round(width/2 + (d/2) + (length * sin(angle2)))
    pendulum.y2 = round(height/4 + (length * cos(angle2)))


def grid():  # Draw a grid behind the pendulum
    for x in range(50, width, 50):
        pygame.draw.lines(background, gray, False, [(x, 0), (x, height)])
        for y in range(50, height, 50):
            pygame.draw.lines(background, gray, False, [(0, y), (width, y)])
    # pygame.draw.circle(background, black, (int(width/2), 50), 5)


def redraw():  # Clean up the screen and start a new grid and new frame of pendulum with new coordinates
    background.fill(white)
    pendulum.draw(background)
    pygame.display.update()


def G_adim(y, t):
    theta1, dtheta1, theta2, dtheta2 = y

    d_theta1 = dtheta1

    dd_theta1 = A*sin(theta1) - B*sin(theta1 - theta2) * \
        (1 - (d / sqrt(d**2 + 2*l**2*(1 - cos(theta1 + theta2))))) - x*dtheta1
    
    d_theta2 = dtheta2

    dd_theta2 = A*sin(theta2) - C*sin(theta2 - theta1) * \
        (1 - (d / sqrt(d**2 + 2*l**2*(1 - cos(theta1 + theta2))))) - x*dtheta2
    
    return np.array([d_theta1, dd_theta1, d_theta2, dd_theta2])


# I start the class with some random coordinates
pendulum = ball((-d/2 + width/2, l + height/4),
                (d/2 + width/2, l + height/4), 5)

# Time points to solve the ODE for
t = np.linspace(0, 100, 1000)

# Solve the ODE
y0 = [theta1_0, dot_theta1_0, theta2_0, dot_theta1_0]
sol = odeint(G_adim, y0, t)

# Extract dd_theta1 and dd_theta2
soltheta1 = sol[:, 0]
sold_theta1 = sol[:, 1]
soltheta2 = sol[:, 2]
sold_theta2 = sol[:, 3]

# Initialize plot
plt.ion()
fig, ax = plt.subplots()
line1, = ax.plot([], [], label='theta1')
line2, = ax.plot([], [], label='theta2')
ax.set_xlabel('Temps')
ax.set_ylabel('Theta')
ax.set_title('Oscillations du système en fonction du temps')
ax.set_xlim(0, 100)
ax.legend()
ax.grid(True)
plt.plot(t, soltheta1, label='theta1', color='blue')
plt.plot(t, soltheta2, label='theta2', color='red')


# Collect data
theta1_data = soltheta1  # Oscillation data for pendulum 1
theta2_data = soltheta2  # Oscillation data for pendulum 2

# Sampling frequency
fs = 1 / (t[1] - t[0])  # Sample rate, inverse of time step

# Perform FFT
theta1_fft = fft(theta1_data)
theta2_fft = fft(theta2_data)

# Frequency axis
freqs = fftfreq(len(theta1_data), 1 / fs)

# Plot FFT results
plt.figure()
plt.plot(freqs, np.abs(theta1_fft), label='Pendulum 1')
plt.plot(freqs, np.abs(theta2_fft), label='Pendulum 2')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude')
plt.title('Frequency Analysis')
plt.xlim(-1,1)  # Zooming on the x-axis, adjust the range as needed
plt.legend()
plt.grid(True)
plt.show()

# Update plot in a loop
t = 0

while not Out:
    print(t)
    for event in pygame.event.get():
            if event.type == pygame.QUIT:
                Out = True  # Si la fenêtre est fermée, sortir de la boucle
            if event.type == pygame.MOUSEBUTTONDOWN:  # Read if you want go out
                pendulum = ball((width/2 - (d/2), height/4 + l),
                                (width/2 + (d/2), height/4 + l), 5)
                acceleration = True

    if acceleration:   # Increase acceleration and damping in the pendulum moviment
        # Increment time
        t += 0.1

        # Solve the ODE
        y0 = [theta1_0, dot_theta1_0, theta2_0, dot_theta2_0]
        sol, info = odeint(G_adim, y0, [0, t], full_output=True)

        # Accédez aux informations détaillées
        print(info)

        # Extract theta1 and theta2
        theta1 = sol[-1, 0]
        theta2 = sol[-1, 2]

        get_path(theta1 + pi, theta2 + pi, l)

        # Update plot data
        line1.set_xdata(np.append(line1.get_xdata(), t))
        line1.set_ydata(np.append(line1.get_ydata(), theta1))
        line2.set_xdata(np.append(line2.get_xdata(), t))
        line2.set_ydata(np.append(line2.get_ydata(), theta2))

        # Redraw the plot
        ax.relim()
        #ax.autoscale_view()

    redraw()

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
