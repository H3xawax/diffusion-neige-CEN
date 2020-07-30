import numpy as np
import matplotlib.pyplot as plt
import progressbar as pb
import matplotlib.colors as mcolors
import time
start_time = time.time()

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


def Hconv(T,epsi, L, Tf,C):
    if T<=(Tf-epsi): return C*T
    if T>(Tf+epsi): return C*T+L
    else: return C*T+L*(T+epsi)/(2*epsi)

def Tconv(H,epsi,L,C):
    if H<=(-C*epsi) : return H/C
    if H>(C*epsi+L) : return (H-L)/C
    else : return epsi*(H-L/2)/(C*epsi+L/2)


def fonte(T,Tf,n,L,C,epsi,rho):
    Ts=T
    solid=1-phase (T, Tf, epsi)
    for i in range(Nx):
        if T[n,i]>=Tf and T[n-1,i]<Tf:
            Ts[n,i]=(C*T[n,i]-min(C*rho*(T[n,i]-Tf),solid[n,i]*L))/C
    return Ts

L=333550. #chaleur latente fiusion de la glace
rho= 917. #masse volumique
K=2.22
C=2060

###################################################
dx=.005 #  metre   /!\ ca ne marche pas avec tous les dx si trop grand on devient absurde
Totprofond=1
Nx=int(Totprofond/dx)
dx2=dx*dx
Tottime=5000
dt=dx2*.45*rho*C/K
Nt= int(Tottime/dt)#nb de pas de temps

###################################################

epsi=.1
Tf=.0
convergence=.001
R=1./(1.+2.*C*epsi/L)

#####################################
bordhaut=20.
Tini=-2.
bordbas=-20.
#################################


T=np.ones((Nt,Nx))*Tini

T[:,0]=bordhaut #voir linspace
T[:,Nx-1]=bordbas

print('############### CFL: ',dt*K/(dx2*rho*C),"###############")
print('Résolution en cours')

A=dt*K/(dx2*rho*C)

Amatrice=np.zeros((Nx,Nx))
Tn_1=np.zeros(Nx)

for i in range(1,Nx-1):
    Amatrice[i, i - 1] = -A
    Amatrice[i, i + 1] = -A
    Amatrice[i, i] = 1 + 2 * A
Amatrice[0, 0] = Amatrice[Nx-1, Nx-1] = 1
b = np.zeros(Nx)

for n in pb.progressbar(range(1,Nt)):
    T[n,:]=np.linalg.solve(Amatrice,T[n-1,:])
    #T=fonte(T, Tf, n, L, C,epsi,rho)

print("--- %s seconds ---" % (time.time() - start_time))

Phase=phase(T,Tf,epsi)
plt.ylabel('Profondeur (m)')
plt.xlabel('Pas de Temps (s)')
xtent = [dt*0 , Nt,  dx*0, Nx]
#print(extent)
norm = mcolors.TwoSlopeNorm(vmin=T.min(), vmax = T.max(), vcenter=0) #pour fixer le 0 au blanc
im=plt.imshow(np.transpose(T),cmap=plt.cm.seismic, norm=norm ,aspect='auto',interpolation='None')
plt.title('CROCUS CFL: '+str(round((dt*K)/(dx2*C),5))+'\n Th: '+ str(bordhaut)+ ' Tb: '+str( bordbas)+ ' Ti: '+str(Tini)+' dt: '+str(round(dt,5))+ " dx: "+str(round(dx,5))+"\n Execution time: "+str(round(time.time() - start_time))+"s")
cb=plt.colorbar()
cb.ax.set_ylabel('Temperature °C', rotation=270)
plt.contourf(np.ma.masked_where(np.transpose(Phase)==0,np.transpose(Phase)),0,hatches=[ '////'], alpha=0)
plt.show()
plt.imshow(np.transpose(Phase),cmap='viridis',aspect='auto',interpolation='None')
plt.contourf(np.ma.masked_where(np.transpose(Phase)==0,np.transpose(Phase)),0,hatches=[ '////'], alpha=0)
plt.colorbar()
plt.title('/!\CROCUS CFL Phase')
plt.show()