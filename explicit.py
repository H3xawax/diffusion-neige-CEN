import numpy as np
import matplotlib.pyplot as plt
import progressbar as pb

def Hconv(T,epsi, L, Tf,C):
    if T<(Tf-epsi): return C*T
    if T>(Tf+epsi): return C*T+L
    else: return C*T+L*(T+epsi)/(2*epsi)

def Tconv(H,epsi,L,C):
    if H<(-C*epsi) : return H/C
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








 #profondeur y a un pb avec le arrange
L=333550. #chaleur latente fiusion de la glace
rho= 917. #masse volumique
Tottime=50000
dt=1  # seconde
Nt= int(Tottime/dt)#nb de pas de temps
Totprofond=1
dx=.01 #  metre   /!\ ca ne marche pas avec tous les dx si trop grand on devient absurde
Nx=int(Totprofond/dx)
dx2=dx*dx
K=2.22
C=2060
epsi=.0001
Tf=0.
convergence=.00000000001
R=1./(1.+2.*C*epsi/L)

bordhaut=20.
bordbas=20.
T=np.ones((Nt,Nx))*-2.
H=np.empty_like(T)
#T=np.ones((Nt,Nx))*np.linspace(bordhaut-5,bordbas,Nx)

T[:,0]=bordhaut #voir linspace
T[:,Nx-1]=bordbas
for t in pb.progressbar(range(Nt)):
    for x in range(Nx):
        H[t,x] = Hconv(T[t,x],epsi, L, Tf,C)
kmax=0
print('############### CFL: ',(dt*K)/(dx2*C),"###############")

for t in pb.progressbar(range(1,Nt)):
    #print('######t: ', t)
    H[t,:]=explicit(H[t-1,:],rho,dx2,C,dt,K,Nx)

for t in pb.progressbar(range(Nt)):
    for x in range(Nx):
        T[t,x] = Tconv(H[t,x],epsi, L,C)
print(':', T[Nt-1,:])

plt.ylabel('Profondeur (m)')
plt.xlabel('Temps (s)')
extent = [dt*0 , Nt,  dx*0, Nx]
#print(extent)
im=plt.imshow(np.transpose(T),cmap='viridis',aspect='auto',interpolation='none')
clb=plt.colorbar(im)
clb.set_label('Temperature')
plt.show()