import numpy as np
import matplotlib.pyplot as plt
import progressbar as pb
import matplotlib.colors as mcolors
import time
start_time = time.time()

def alpha (T, Tf, epsi,L):
    if T<Tf-epsi: return 0
    if T>Tf+epsi: return -L
    else : return -L/2

def beta (T, Tf, epsi,L,C,lambd):
    if T<Tf-epsi: return C+lambd
    if T>Tf+epsi: return C+lambd
    else : return C+lambd+L/(2*epsi)

def Phi_1 (T, Tf, epsi,L,C):
    if T<Tf-epsi: return C*T
    if T>Tf+epsi: return C*T+L
    else : return C*T+L*(T+epsi)/(2*epsi)


 #profondeur y a un pb avec le arrange
L=333550. #chaleur latente fiusion de la glace
rho= 917. #masse volumique
K=2.22
C=2060

###############################################
# Tottime=500000
# dt=10000# seconde
# Nt= int(Tottime/dt)#nb de pas de temps
# Totprofond=1
# dx=.01#np.sqrt(dt*K/(C*.43)) #  metre   /!\ ca ne marche pas avec tous les dx si trop grand on devient absurde
# Nx=int(Totprofond/dx)
# dx2=dx*dx
# print(dt)
# print(dx2)
# print(Nt)
# print(Nx)


###################################################
dx=.025 #  metre   /!\ ca ne marche pas avec tous les dx si trop grand on devient absurde
Totprofond=1
Nx=int(Totprofond/dx)
dx2=dx*dx
Tottime=50000
dt=dx2*C*.45/K
Nt= int(Tottime/dt)#nb de pas de temps

###################################################

epsi=.5
Tf=.0
lambd = 2*dt*K/(dx2*rho)
#print(dx)
#print((1+lambd/C))
#print(L/(C*(1+lambd/C)))

#print((L/C)/(2+(2*lambd/C)))

#####################################
bordhaut=20.
Tini=-2.
bordbas=-20.
#################################


T=np.ones((Nt,Nx))*Tini
#T=np.ones((Nt,Nx))*np.linspace(bordhaut-5,bordbas,Nx)

T[:,0]=bordhaut #voir linspace
T[:,Nx-1]=bordbas
print('############### CFL: ',(dt*K)/(dx2*C),"###############")

for t in pb.progressbar(range(Nt-1)):
    for i in range(1,Nx-1):
        T[t+1,i]=( Phi_1(T[t,i],Tf, epsi,L,C) + (dt*K)*(T[t,i-1]+T[t,i+1])/(dx2*rho) + alpha(T[t,i], Tf, epsi,L))/beta(T[t,i],Tf, epsi,L,C,lambd)

print("--- %s seconds ---" % (time.time() - start_time))

plt.ylabel('Profondeur (m)')
plt.xlabel('Temps (s)')
extent = [dt*0 , Nt,  dx*0, Nx]
#print(extent)
norm = mcolors.TwoSlopeNorm(vmin=T.min(), vmax = T.max(), vcenter=0) #pour fixer le 0 au blanc
im=plt.imshow(np.transpose(T),cmap=plt.cm.seismic, norm=norm ,aspect='auto',interpolation='none')
plt.title('/!\VOLLER2 CFL: '+str(round((dt*K)/(dx2*C),5))+'\n Th: '+ str(bordhaut)+ ' Tb: '+str( bordbas)+ ' Ti: '+str(Tini)+' dt: '+str(round(dt,5))+ " dx: "+str(round(dx,5))+"\n Execution time: "+str(round(time.time() - start_time))+"s")
plt.colorbar()
#cax.set_label('Temperature')
plt.show()