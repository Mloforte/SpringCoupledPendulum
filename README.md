# SpringCoupledPendulum

This is a simple simulation of a spring coupled pendulum. the physics of the system is described by the following equations:

![equation](https://latex.codecogs.com/svg.image?\[\begin{align*}d_{\theta_1}&=\dot{\theta}_1\\\ddot{\theta}_1&=A\sin(\theta_1)-B\sin(\theta_1-\theta_2)\left(1-\left(\frac{d}{\sqrt{d^2&plus;2l^2(1-\cos(\theta_1&plus;\theta_2))}}\right)\right)-x\dot{\theta}_1\\\end{align*}\])

![equation](https://latex.codecogs.com/svg.image?\[\begin{align*}d_{\theta_2}&=\dot{\theta}_2\\\ddot{\theta}_2&=A\sin(\theta_2)-C\sin(\theta_2-\theta_1)\left(1-\left(\frac{d}{\sqrt{d^2&plus;2l^2(1-\cos(\theta_1&plus;\theta_2))}}\right)\right)-x\dot{\theta}_2\\\end{align*}\])

![alt text](Image/image.png)

The spring is connected to the two masses and the spring constant is k. the length of the pendulum is l and the masses are m1 and m2. the angle of the pendulum is theta1 and theta2. the angular velocity of the pendulum is dot_theta1 and dot_theta2. the distance between the two masses is d. the acceleration due to gravity is g. The damping coefficient is x to account for air resistance.

The two masses are connected to a pivot point by a massless rod. The rod is assumed to be rigid and the masses are assumed to be point masses. The rod is assumed to be massless and the spring is assumed to be massless.  

A problem occured after a certain amount of simulation time, I think it's because of the equations of motion.It says :'Excess work done on this call (perhaps wrong Dfun type).' line 153[sol, info = odeint(G_adim, y0, [0, t], full_output=True)] . I would appreciate it if someone could help me with this.