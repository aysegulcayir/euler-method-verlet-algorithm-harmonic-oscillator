# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 01:24:22 2020

@author: pc
"""

from numpy import zeros, linspace, cos, array
import matplotlib.pyplot as plt

k=5
m=2
dt =0.1 
T = 20
N_t = int(round(T/dt))
t = linspace(0, N_t*dt, N_t+1)


u = zeros(N_t+1)
v = zeros(N_t+1)
ek = zeros(N_t+1)
ep = zeros(N_t+1)


# Initial condition
X_0 = 0.01
u[0] = X_0
v[0] = 0
ek[0] = 0
ep[0] = 0

print("\n analytical solutions:\n")
print(X_0*cos(((5/2)**(0.5))*t))

# Step equations forward in time
print("\n numerical solutions with forward euler:\n")
for n in range(N_t):
    v[n+1] = v[n] - dt*(k/m)*u[n]
    u[n+1] = u[n] + dt*v[n+1]
    v[n+1] = v[n] - dt*(k/m)*u[n+1]
    print( u[n+1] , v[n+1])

    ek[n] = (0.5) * m * ((v[n]) ** 2)
    ep[n] = (0.5) * k * ((u[n]) ** 2)

fig = plt.figure()
l1, l2 = plt.plot(t, u, 'b-', t, X_0*cos(((5/2)**(0.5))*t), 'r--')
fig.legend((l1, l2), ('forwardEuler', 'exact'), 'upper left')
plt.xlabel('t')
plt.show()

fig = plt.figure()
lk,lp = plt.plot(t, ek, t, ep, 'r--')
fig.legend((lk,lp), ('kinetic_energy','potential_energy'), 'upper left')
plt.xlabel('t')
plt.show()

print("\n numerical solutions with verlet:\n")

for n in range(N_t):
    
    u[n+1] = u[n] + dt*v[n]
    v[n+1] = v[n] - dt*(k/m)*u[n+1]
    print( u[n+1] , v[n+1])

fig1 = plt.figure()
l11, l22 = plt.plot(t, u, 'b-', t, X_0*cos(((5/2)**(0.5))*t), 'r--')
fig1.legend((l11, l22), ('verlet', 'exact'), 'upper left')
plt.xlabel('t')
plt.show()


plt.savefig('tmp.pdf'); plt.savefig('tmp.png')