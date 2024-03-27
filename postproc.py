 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 16:07:20 2020

@author: louisalsteens
"""

import numpy as np
import matplotlib.pyplot as plt
import readField as rf
import csv

import readField as rf
import readIni

plt.close('all')



#%%

simu = 'CBL_01_avg'
time = 50000
sim = 'CBL_01'

xT,T = rf.readScalar(sim, 'theta0', time)
xu,u = rf.readVel(sim, time,comp=0);
xv,v = rf.readVel(sim, time,comp=1);
xw,w = rf.readVel(sim, time,comp=2);
#xv,v = rf.readStats(simu,comp=1);



u_grid = [[[0 for _ in range(len(u[0][0]))] for _ in range(len(u[0]))] for _ in range(len(u))]

# INTERPOLATION DE U 
for i in range(len(u_grid)) : 
    for j in range(len(u_grid[0])):
        for k in range(len(u_grid[0][0])): 
            
            
            if i == 0 and k == (len(u_grid[0][0]) - 1 ):
                u_grid[i][j][k] = u[i][j][k]
                    
                
            if i == 0 and k != (len(u_grid[0][0]) - 1 ):
                u_grid[i][j][k] = (u[i][j][k+1] + u[i][j][k]) / 2
                    
                    
            if k == (len(u_grid[0][0]) - 1) and i != 0 : 
                u_grid[i][j][k] = (u[i-1][j][k] + u[i][j][k]) / 2 
            
            
            if k != (len(u_grid[0][0]) - 1) and i != 0 :
                u_grid[i][j][k] = (u[i-1][j][k+1] + u[i-1][j][k] + u[i][j][k+1] + u[i][j][k]) / 4 
      
u_grid_array = np.array(u_grid)

"""
v_grid = [[[0 for _ in range(len(u[0][0]))] for _ in range(len(u[0]))] for _ in range(len(u))]
w_grid = [[[0 for _ in range(len(u[0][0]))] for _ in range(len(u[0]))] for _ in range(len(u))]

# INTERPOLATION DE V 
for i in range(len(u_grid)) : 
    for j in range(len(u_grid[0])):
        for k in range(len(u_grid[0][0])): 
            
            
            if i == 0 and j == 0 :
                v_grid[i][j][k] = v[i][j][k]
                
                
            if i == 0 and j != 0 :
                v_grid[i][j][k] = (v[i][j-1][k] + v[i][j][k]) / 2
            
            
            if i != 0 and j == 0 : 
                v_grid[i][j][k] = (v[i-1][j][k] + v[i][j][k]) / 2
            
            
            if i != 0 and j != 0 :
                v_grid[i][j][k] = (v[i-1][j-1][k] + v[i-1][j][k] + v[i][j-1][k] + v[i][j][k]) / 4


 # INTERPOLATION DE W
for i in range(len(u_grid)) : 
    for j in range(len(u_grid[0])):
        for k in range(len(u_grid[0][0])): 
            
            
            if j == 0 and k == (len(u_grid[0][0]) - 1 ) : 
                w_grid[i][j][k] = w[i][j][k] 
            
            
            if j != 0 and k == (len(u_grid[0][0]) - 1 ) : 
                w_grid[i][j][k] = (w[i][j-1][k] + w[i][j][k]) / 2
                
                
            if j == 0 and k != (len(u_grid[0][0]) - 1 ) : 
                w_grid[i][j][k] = (w[i][j-1][k] + w[i][j][k+1]) / 2
            
            
            if j != 0 and k != (len(u_grid[0][0]) - 1 ) : 
                w_grid[i][j][k] = (w[i][j-1][k] +  w[i][j-1][k+1] + w[i][j][k+1] + w[i][j][k]) / 4
            
          
"""
            
"""           
u_aplati = u.flatten()
v_aplati = v.flatten()
w_aplati = w.flatten()


np.savetxt('u_aplati.csv', u_aplati, delimiter=',')
np.savetxt('v_aplati.csv', v_aplati, delimiter=',')
np.savetxt('w_aplati.csv', w_aplati, delimiter=',')
"""


#aplati1 = xu[0].flatten()
#aplati2 = xu[1].flatten()
#aplati3 = xu[2].flatten()

# Concaténer les tableaux aplatis
#tableau_final = np.concatenate([aplati1, aplati2, aplati3])

# Sauvegarder le tableau final dans un fichier CSV
#np.savetxt('grid.csv', tableau_final, delimiter=',')





#%%
 

f=0;
#%%
#Plotting the results
"""
#Umean velocity
f=f+1;
level = np.linspace(300.0,300.7,100)
plt.figure(num=f,figsize=(6,5));
plt.contourf(xT[1][3:-3],xT[2][3:-3],T[50,3:-3,3:-3], levels=level, cmap=plt.cm.seismic);
plt.colorbar();
plt.legend(loc="upper right")
ax = plt.gca()
ax.set_xlim(0, 8000)
ax.set_ylim(0, 8000)
# Set axis limits and equals
#ax = plt.gca()
#ax.set(xlim=(-2,6),ylim=(-2,2))
#ax.set_aspect('equal', adjustable='box', anchor='C') #equivalent to plt.axis('scaled')

# labels
plt.xlabel('Temperature field')
plt.ylabel('Y')
plt.title('Temperature field at z/zi=0.3')
plt.grid(color='k', linestyle=':', linewidth=0.5)
"""
#Umean velocity
data14 = np.loadtxt('/Users/Simon/Documents/Aaa_Thesis/git_jsbsim/jsbsim/Zzz_North.txt')
for i in range(len(data14)):
    data14[i] += 4000
    
data24 = np.loadtxt('/Users/Simon/Documents/Aaa_Thesis/git_jsbsim/jsbsim/Zzz_East.txt')
for i in range(len(data24)):
    data24[i] += 1000


f=f+1;
level = np.linspace(-4,4,100)
plt.figure(num=f,figsize=(6,5));

#BIEN FAIRE ATTENTION A LA HAUTEUR 97 QUELLE CORRESPONDE BIEN A LA HAUTEUR QU ON VEUT REGARDER 

plt.contourf(xT[1][3:-3],xT[2][3:-3],u_grid_array[97,3:-3,3:-3], levels=level, cmap=plt.cm.seismic); 
plt.colorbar();

#RAJOUT DU CHEMIN 
plt.plot(data24, data14,color = 'black', label = 'Trajectoire Planneur', linewidth=3.0)
plt.scatter(data24[0], data14[0], color='green', label='Début', s=70, zorder=5)
plt.scatter(data24[-1], data14[-1], color='red', label='Fin', s=70, zorder=5)

plt.legend(loc="upper right")
ax = plt.gca()
ax.set_xlim(0, 4000)
ax.set_ylim(0, 8000)
# Set axis limits and equals
#ax = plt.gca()
#ax.set(xlim=(-2,6),ylim=(-2,2))
#ax.set_aspect('equal', adjustable='box', anchor='C') #equivalent to plt.axis('scaled')

# labels
plt.xlabel('U Velocity')
plt.ylabel('Y')
plt.title('Mean velocity profile U')
plt.grid(color='k', linestyle=':', linewidth=0.5)

# labels
plt.xlabel('Temperature field')
plt.ylabel('Y')
plt.title('Temperature field at z/zi=0.3')
plt.grid(color='k', linestyle=':', linewidth=0.5)

#%%
#Umean velocity
"""
U = u[:,90,3:-3]
V = v[:,90,3:-3]

f=f+1;
level = np.linspace(300.0,300.8,100)
plt.figure(num=f,figsize=(16,5));
plt.contourf(xT[2][3:-3],xT[0][:],(T[:,90,3:-3]), levels=level, cmap=plt.cm.seismic);
plt.colorbar();
#plt.contourf(xT[1][3:-3],xT[0][3:-3],v[:,:,90], cmap=plt.cm.seismic);
#plt.quiver(Y[::10,::10],X[::10,::10],(V[::10,::10]).T,(U[::10,::10]).T,color='white')
plt.legend(loc="upper right")
ax = plt.gca()
ax.set_xlim(0, 8000)
ax.set_ylim(0,2400)
# Set axis limits and equals
#ax = plt.gca()
#ax.set(xlim=(-2,6),ylim=(-2,2))
#ax.set_aspect('equal', adjustable='box', anchor='C') #equivalent to plt.axis('scaled')
# labels
plt.xlabel('U Velocity')
plt.ylabel('Y')
plt.title('Mean velocity profile U')
plt.grid(color='k', linestyle=':', linewidth=0.5)


# labels
plt.xlabel('Temperature field')
plt.ylabel('Y')
plt.title('Temperature field at z/zi=0.3')
plt.grid(color='k', linestyle=':', linewidth=0.5)


"""
#Umean velocity
data34 = np.loadtxt('/Users/Simon/Documents/Aaa_Thesis/git_jsbsim/jsbsim/Zzz_Alt.txt')


f=f+1;
level = np.linspace(-4,4,100)
plt.figure(num=f,figsize=(16,5));

# ICI LARGEUR 50 : VOIR AVEC PRINT(XU[2][50])

plt.contourf(xT[2][3:-3],xT[0][:],(u_grid_array[:,3:-3,50]), levels=level, cmap=plt.cm.seismic);
#plt.contourf(xT[2][3:-3],xT[0][:],(u[:,3:-3,50]), levels=level, cmap=plt.cm.seismic);
print(xu[2][50])
plt.colorbar();

plt.plot(data14, data34,color = 'black', label = 'Trajectoire Planneur', linewidth=3.0)
plt.scatter(data14[0], data34[0], color='green', label='Début', s=70, zorder=5)
plt.scatter(data14[-1], data34[-1], color='red', label='Fin', s=70, zorder=5)

#plt.contourf(xT[1][3:-3],xT[0][3:-3],v[:,:,90], cmap=plt.cm.seismic);
#plt.quiver(Y[::10,::10],X[::10,::10],(v[::10,::10,90]).T,(u[::10,::10,90]).T,color='white')
plt.legend(loc="upper right")
ax = plt.gca()
ax.set_xlim(0, 8000)
ax.set_ylim(0,2400)
plt.axis('equal')

# labels
plt.xlabel('U Velocity')
plt.ylabel('Y')
plt.title('Mean velocity profile U first interp')
plt.grid(color='k', linestyle=':', linewidth=0.5)

