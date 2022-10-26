#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 12:54:51 2022

@author: elenamunozrivas
"""

import matplotlib.pyplot as plt
import numpy as np

#Movement of the first Paricle in the Penning Trap:
PeningTrap1 = np.loadtxt("PenningTrap1.txt")

#Import the values of time and x in Lists:
time = PeningTrap1[:,0]
z = PeningTrap1[:,1]

#Plot a motion of x direction as a function of time
plt.plot(time,z, color = "darkgreen")
plt.scatter(time[0],z[0], color = "darkgreen")
plt.scatter(time[len(time)-1],z[len(z)-1], color = "darkgreen")
plt.xlabel("t (μs)")
plt.ylabel("z (μm)")
plt.gca().set_aspect('equal')
plt.grid(True)
plt.legend()
plt.savefig('time,z_PT1.pdf')

plt.show()

#=============================================================================

#Movement of the two particles with interactions:
PeningTrap2 = np.loadtxt("PenningTrap2.txt")

#Import the values of Penning Trap 2 in Lists:
x1 = PeningTrap2[:,1]
y1 = PeningTrap2[:,2]
z1 = PeningTrap2[:,3]
vx1 = PeningTrap2[:,4]
vz1 = PeningTrap2[:,6]
x2 = PeningTrap2[:,7]
y2 = PeningTrap2[:,8]
z2 = PeningTrap2[:,9]
vx2 = PeningTrap2[:,10]
vz2 = PeningTrap2[:,12]

#Plot the values of y axis against x direction:
plt.plot(x1,y1, label="Particle 1", color = "darkgreen")
plt.scatter(x1[0],y1[0], color = "darkgreen")
plt.scatter(x1[len(x1)-1],y1[len(y1)-1], color = "darkgreen")
plt.plot(x2,y2, label = "Particle 2", color = "darkorange")
plt.scatter(x2[0],y2[0], color = "darkorange")
plt.scatter(x2[len(x2)-1],y2[len(y2)-1], color = "darkorange")
plt.xlabel("x (μm)")
plt.ylabel("y (μm)")
plt.gca().set_aspect('equal')
plt.grid(True)
plt.legend()
plt.savefig('x,y_with.pdf')

plt.show()

#=============================================================================


#Plot a motion of v_x as a function of x for the two particles:
plt.plot(x1,vx1, label="Particle 1", color = "darkgreen")
plt.scatter(x1[0],vx1[0], color = "darkgreen")
plt.scatter(x1[len(x1)-1], vx1[len(vx1)-1], color = "darkgreen")
plt.plot(x2,vx2, label = "Particle 2", color = "darkorange")
plt.scatter(x2[0],vx2[0], color = "darkorange")
plt.scatter(x2[len(x2)-1],vx2[len(vx2)-1], color = "darkorange")
plt.xlabel("x (μm)")
plt.ylabel("$v_x$ (μm/μs)")
plt.gca().set_aspect('equal')
plt.grid(True)
plt.legend(loc='upper right')
plt.savefig('x,vx_with.pdf')

plt.show()

#=============================================================================


#Plot a motion of v_z as a function of z for the two particles:
plt.plot(z1,vz1, label="Particle 1", color = "darkgreen")
plt.scatter(z1[0],vz1[0], color = "darkgreen")
plt.scatter(z1[len(z1)-1], vz1[len(vz1)-1], color = "darkgreen")
plt.plot(z2,vz2, label="Particle 2", color = "darkorange")
plt.scatter(z2[0],vz2[0], color = "darkorange")
plt.scatter(z2[len(z2)-1], vz2[len(vz2)-1], color = "darkorange")
plt.xlabel("z (μm)")
plt.ylabel("$v_z$ (μm/μs)")
plt.gca().set_aspect('equal')
plt.grid(True)
plt.legend(loc='upper right')
plt.savefig('z,vz_with.pdf')

plt.show()

#=============================================================================


#3D Plot of the trajectory of two particles in the space
#without interactions:
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111, projection='3d')
ax.plot(x1,y1,z1, label="Particle 1"  ,color = "darkgreen")
ax.plot(x2,y2,z2, label="Particle 2"  ,color = "darkorange")

ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z');
ax.legend()
plt.savefig('3D_axis_with.pdf')

plt.show()

#=============================================================================


#Movement of the two particles with interactions:
PeningTrap2_without = np.loadtxt("PenningTrap2_without.txt")

#Import the values of Penning Trap 2 without iterations in Lists:
x1 = PeningTrap2_without[:,1]
y1 = PeningTrap2_without[:,2]
z1 = PeningTrap2_without[:,3]
vx1 = PeningTrap2_without[:,4]
vz1 = PeningTrap2_without[:,6]
x2 = PeningTrap2_without[:,7]
y2 = PeningTrap2_without[:,8]
z2 = PeningTrap2_without[:,9]
vx2 = PeningTrap2_without[:,10]
vz2 = PeningTrap2_without[:,12]

#Plot a motion of y axis against of x directions for the two particles
#without particle interactions:
plt.plot(x1,y1,  label="Particle 1", color = "darkgreen")
plt.scatter(x1[0],y1[0], color = "darkgreen")
plt.scatter(x1[len(x1)-1], y1[len(y1)-1], color = "darkgreen")

plt.plot(x2,y2,  label="Particle 2", color = "darkorange")
plt.scatter(x2[0],y2[0], color = "darkorange")
plt.scatter(x2[len(x2)-1], y2[len(y2)-1], color = "darkorange")
plt.xlabel("x (μm)")
plt.ylabel("y (μm)")
plt.gca().set_aspect('equal')
plt.grid(True)
plt.legend()
plt.savefig('x,y_without.pdf')

plt.show()

#=============================================================================


#Plot a motion of v_x as a function of x for the two particles
#without particle interactions:
plt.plot(x1,vx1, label="Particle 1"  ,color = "darkgreen")
plt.scatter(x1[0],vx1[0], color = "darkgreen")
plt.scatter(x1[len(x1)-1], vx1[len(vx1)-1], color = "darkgreen")

plt.plot(x2,vx2, label="Particle 2", color = "darkorange")
plt.scatter(x2[0],vx2[0], color = "darkorange")
plt.scatter(x2[len(x2)-1], vx2[len(vx2)-1], color = "darkorange")

plt.xlabel("x (μm)")
plt.ylabel("$v_x$ (μm/μs)")
plt.gca().set_aspect('equal')
plt.grid(True)
plt.legend(loc='upper right')
plt.savefig('x,vx_without.pdf')

plt.show()

#=============================================================================


#Plot a motion of v_z as a function of z direction for the two particles
#without particle interactions:
plt.plot(z1,vz1, label="Particle 1"  ,color = "darkgreen")
plt.scatter(z1[0],vz1[0], color = "darkgreen")
plt.scatter(z1[len(z1)-1], vz1[len(vz1)-1], color = "darkgreen")

plt.plot(z2,vz2, label="Particle 2"  ,color = "darkorange")
plt.scatter(z2[0],vz2[0], color = "darkorange")
plt.scatter(z2[len(z2)-1], vz2[len(vz2)-1], color = "darkorange")

plt.xlabel("z (μm)")
plt.ylabel("$v_z$ (μm/μs)")
plt.gca().set_aspect('equal')
plt.grid(True)
plt.legend(loc='upper right')
plt.savefig('z,vz_without.pdf')

plt.show()

#=============================================================================

#3D Plot of the trajectory of two particles in the space
#without interactions:
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111, projection='3d')
ax.plot(x1,y1,z1, label="Particle 1"  ,color = "darkgreen")
ax.plot(x2,y2,z2, label="Particle 2"  ,color = "darkorange")

ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z');
ax.legend()
plt.savefig('3D_axis_without.pdf')

plt.show()

#=============================================================================

#Load the txt with the relative errors
datos1= np.loadtxt('relative_errorn1_eu_rku.txt')
datos2= np.loadtxt('relative_errorn2_eu_rku.txt')
datos3= np.loadtxt('relative_errorn3_eu_rku.txt')
datos4= np.loadtxt('relative_errorn4_eu_rku.txt')

#Import the values of the relatives errors for each n:
t1 = datos1[:,0]
err1_rk4 = datos1[:,1]
err1_eu = datos1[:,2]

t2= datos2[:,0]
err2_rk4=datos2[:,1]
err2_eu=datos2[:,2]

t3= datos3[:,0]
err3_rk4=datos3[:,1]
err3_eu=datos3[:,2]

t4= datos4[:,0]
err4_rk4=datos4[:,1]
err4_eu=datos4[:,2]

#Plot the relatives errors as functions of time:
plt.plot(t1,err1_rk4, label='n=4000',color='darkorange')
plt.plot(t1,err1_eu,ls='dotted',color='darkorange')
plt.plot(t2,err2_rk4, label='n=8000',color='darkgreen')
plt.plot(t2,err2_eu,color='darkgreen', ls='dotted')
plt.plot(t3,err3_rk4, label='n=16000',color='blue')
plt.plot(t3,err3_eu,color='blue', ls='dotted')
plt.plot(t4,err4_rk4,label='n=32000', color='purple')
plt.plot(t4,err4_eu, color='purple', ls='dotted')

plt.xlabel('t(μs)')
plt.ylabel('$r_{err}$')
plt.grid()
plt.legend()
plt.savefig('relative_error.pdf')

plt.show()

#=============================================================================

#Eliminating the first values ​​so that there is no error when doing the log:
#Do the log of  the relatives errors for each n:
t1 = t1[1:]
err1_log_rk4 = np.log(err1_rk4[1:])
err1_log_eu = np.log(err1_eu[1:])
t2 = t2[1:]
err2_log_rk4=np.log(err2_rk4[1:])
err2_log_eu=np.log(err2_eu[1:])
t3 = t3[1:]
err3_log_rk4=np.log(err3_rk4[1:])
err3_log_eu=np.log(err3_eu[1:])
t4 = t4[1:]
err4_log_rk4=np.log(err4_rk4[1:])
err4_log_eu=np.log(err4_eu[1:])


#Plot the relatives errors as functions of time:
plt.plot(t1,err1_log_rk4, label='n=4000',color='darkorange')
plt.plot(t1,err1_log_eu,ls='dotted',color='darkorange')
plt.plot(t2,err2_log_rk4, label='n=8000',color='darkgreen')
plt.plot(t2,err2_log_eu,color='darkgreen', ls='dotted')
plt.plot(t3,err3_log_rk4, label='n=16000',color='blue')
plt.plot(t3,err3_log_eu,color='blue', ls='dotted')
plt.plot(t4,err4_log_rk4,label='n=32000', color='purple')
plt.plot(t4,err4_log_eu, color='purple', ls='dotted')

plt.xlabel('t(μs)')
plt.ylabel('$r_{err}$')
plt.grid()
plt.legend(loc='lower right', prop={'size': 8.5})
plt.savefig('relative_error_log.pdf')

plt.show()

#=============================================================================

#Load the txt with the fraction of particles values for each frecuencies:
datosf1= np.loadtxt('PenningTrap100f1.txt')
datosf4= np.loadtxt('PenningTrap100f4.txt')
datosf7= np.loadtxt('PenningTrap100f7.txt')

#Import the values in Lists:
w_v1 = datosf1[:,0]
f_part1 = datosf1[:,1]
w_v4 = datosf4[:,0]
f_part4 = datosf4[:,1]
w_v7 = datosf7[:,0]
f_part7 = datosf7[:,1]


#Plot the fraction of particles that are still trapped after 500µs
#as a function of the applied angular frequency w_v:
plt.plot(w_v1,f_part1, label='f=0.1',color='darkorange')
plt.plot(w_v4,f_part4, label='f=0.4',color='darkgreen')
plt.plot(w_v7,f_part7, label='f=0.7',color='blue')


plt.xlabel('$w_v (MHz) $')
plt.ylabel('Fraction of particles')
plt.grid()
plt.legend(loc='lower left', prop={'size': 8.5})
plt.savefig('fractpart.pdf')

plt.show()


#=============================================================================

#Load the txt with the fraction of particles values with zoom with and without
#particle interactions:
datosf_with= np.loadtxt('PenningTrap100zoom1.txt')
datosf_without= np.loadtxt('PenningTrap100zoom2.txt')
datosf7= np.loadtxt('PenningTrap100f7.txt')

#Import the values in Lists:
w_v_with = datosf_with[:,0]
f_part_with = datosf_with[:,1]
w_v_without = datosf_without[:,0]
f_part_without = datosf_without[:,1]


#Plot the fraction of trapped particles versus frequency with and without
#Coulomb interactions :
plt.plot(w_v_with,f_part_with, label='f=with interactions',color='darkorange')
plt.plot(w_v_without,f_part_without, label='f=without interactions',color='darkgreen')


plt.xlabel('$w_v (MHz) $')
plt.ylabel('Fraction of particles')
plt.grid()
plt.legend(loc='lower right', prop={'size': 7.5})
plt.savefig('fractpart_intec.pdf')

plt.show()




