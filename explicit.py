import numpy as np
import matplotlib.pyplot as plt
import progressbar as pb
import matplotlib.colors as mcolors
import time
start_time = time.time()

def Hconv(T,epsi, L, Tf,C):
    if T<=(Tf-epsi): return C*T
    if T>(Tf+epsi): return C*T+L
    else: return C*T+L*(T+epsi)/(2*epsi)

def Tconv(H,epsi,L,C):
    if H<=(-C*epsi) : return H/C
    if H>(C*epsi+L) : return (H-L)/C
    else : return epsi*(H-L/2)/(C*epsi+L/2)


def explicit(Hj,rho,dx2,C,dt,K,Nx):
    Hj1=np.ones(len(Hj))
    Hj1[0]=Hj[0]
    Hj1[-1] = Hj[-1]
    for i in range(1,Nx-1):
        Hj1[i]=Hj[i]+(Tconv(Hj[i-1],epsi,L,C)-2*Tconv(Hj[i],epsi,L,C)+Tconv(Hj[i+1],epsi,L,C))*(dt*K)/(rho*dx2)
        #print(i,'tada:',(dt*K)/(rho*dx2*C))
    return Hj1

def phase (T, Tf, epsi):
    Phase=np.empty_like(T,dtype=float)
    for i in range(np.shape(T)[0]):
        for j in range(np.shape(T)[1]):
            if T[i,j]<=(Tf-epsi):
                Phase[i,j]=0.
            if T[i,j]>=(Tf+epsi):
                Phase[i,j]=1.
            if (T[i,j]<(Tf+epsi) and T[i,j]>(Tf-epsi)):
                Phase[i,j]=(T[i,j]+Tf+epsi)/(2*epsi)
    return Phase

L=333550. #chaleur latente fiusion de la glace
rho= 917. #masse volumique
K=2.22
C=2060

###################################################
dx=.025 #  metre   /!\ ca ne marche pas avec tous les dx si trop grand on devient absurde
Totprofond=1
Nx=int(Totprofond/dx)
dx2=dx*dx
Tottime=50000
dt=dx2*C*45/K
Nt= int(Tottime/dt)#nb de pas de temps

###################################################

epsi=1
Tf=.0
convergence=.001
R=1./(1.+2.*C*epsi/L)

#####################################
bordhaut=20.
Tini=-2.
bordbas=-20.
#################################


T=np.ones((Nt,Nx))*Tini
H=np.empty_like(T)
#T=np.ones((Nt,Nx))*np.linspace(bordhaut-5,bordbas,Nx)

T[:,0]=bordhaut #voir linspace
T[:,Nx-1]=bordbas

print('Transformation de la matrice en Entalpie')
for t in pb.progressbar(range(Nt)):
    for x in range(Nx):
        H[t,x] = Hconv(T[t,x],epsi, L, Tf,C)
kmax=0
print()
print('############### CFL: ',(dt*K)/(dx2*C),"###############")
print('Résolution de la diffusion en Entalpie')

for t in pb.progressbar(range(1,Nt)):
    #print('######t: ', t)
    H[t,:]=explicit(H[t-1,:],rho,dx2,C,dt,K,Nx)

print('Transformation de la matrice en Temperature')
for t in pb.progressbar(range(Nt)):
    for x in range(Nx):
        T[t,x] = Tconv(H[t,x],epsi, L,C)
#print(':', T[Nt-1,:])

print("--- %s seconds ---" % (time.time() - start_time))
Phase=phase(T,Tf,epsi)

plt.ylabel('Profondeur (m)')
plt.xlabel('Temps (s)')
xtent = [dt*0 , Nt,  dx*0, Nx]
#print(extent)
norm = mcolors.TwoSlopeNorm(vmin=T.min(), vmax = T.max(), vcenter=0) #pour fixer le 0 au blanc
im=plt.imshow(np.transpose(T),cmap=plt.cm.seismic, norm=norm ,aspect='auto',interpolation='None')
plt.title('EXPLICIT CFL: '+str(round((dt*K)/(dx2*C),5))+'\n Th: '+ str(bordhaut)+ ' Tb: '+str( bordbas)+ ' Ti: '+str(Tini)+' dt: '+str(round(dt,5))+ " dx: "+str(round(dx,5))+"\n Execution time: "+str(round(time.time() - start_time))+"s")
cb=plt.colorbar()
cb.ax.set_ylabel('Temperature °C', rotation=270)
plt.contourf(np.transpose(Phase),1,hatches=['', '////'], alpha=0,aspect='auto',interpolation='None')
plt.show()
