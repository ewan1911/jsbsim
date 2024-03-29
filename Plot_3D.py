# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 17:15:29 2024

@author: test
"""

from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt


data1 = np.loadtxt('/Users/Simon/Documents/Aaa_Thesis/git_jsbsim/jsbsim/Zzz_North.txt')
data2 = np.loadtxt('/Users/Simon/Documents/Aaa_Thesis/git_jsbsim/jsbsim/Zzz_East.txt')
data3 = np.loadtxt('/Users/Simon/Documents/Aaa_Thesis/git_jsbsim/jsbsim/Zzz_Alt.txt')

data_aR = np.loadtxt('/Users/Simon/Documents/Aaa_Thesis/git_jsbsim/jsbsim/Zzz_aR.txt')
data_T = np.loadtxt('/Users/Simon/Documents/Aaa_Thesis/git_jsbsim/jsbsim/Zzz_Time.txt')
dara_roll = np.loadtxt('/Users/Simon/Documents/Aaa_Thesis/git_jsbsim/jsbsim/Zzz_Roll.txt')

ax = plt.axes(projection='3d')


ax.set_xlabel('East [m]')
ax.set_ylabel('North [m]')
ax.set_zlabel('Altitude [m]')


ax.set_title('3D : Chemin du planneur dans le ciel ')


ax.plot3D(data2, data1, data3)
ax.scatter3D(data2[0], data1[0], data3[0], color = "green", label ="start", s=70)
ax.scatter3D(data2[-1], data1[-1], data3[-1], color = "red",label ="finish",s=70)

plt.legend()
plt.show()
#--------------------------------------------#


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))


ax1.plot(data2, data1)
ax1.grid(True)
ax1.set_xlabel('East [m]')
ax1.set_ylabel('North [m]')
ax1.set_title(' 2D : North en fonction de East')
#ax1.scatter(data2[0], data1[0], color='green', label='Début', s=70, zorder=5)
#ax1.scatter(data2[-1], data1[-1], color='red', label='Fin', s=70, zorder=5)
ax1.scatter(1000.0, 400.0, color = 'purple', label = 'target', s = 70)


ax2.plot(data1, data3)
ax2.grid(True)
ax2.set_xlabel('North [m]')
ax2.set_ylabel('Altitude [m]')
ax2.set_title(' 2D : Altitude en fonction de North')
#ax2.scatter(data1[0], data3[0], color='green', label='Début', s=70, zorder=5)
#ax2.scatter(data1[-1], data3[-1], color='red', label='Fin', s=70, zorder=5)


plt.tight_layout()

ax1.axis('equal')
ax2.axis('equal')

plt.show()

plt.figure("aileron R position overtime")

plt.plot(data_T, data_aR)
plt.plot(data_T, dara_roll)
plt.plot()
plt.grid()
plt.show()