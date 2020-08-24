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


def fonte(T,Tf,n,L,C,PhaseSL,rho,dx,Nx):
    T=T+273.15
    Ts=T
    Tf=Tf+273.15
    for i in range(1,Nx):
        if T[n,i]>Tf and PhaseSL[i]!=0:
            #print('Phase',PhaseSL[i])
            Etot = C * T[n, i] * rho * dx
            #print('Etot',Etot)
            Efonte = min((C * (T[n, i] - Tf) * rho *dx),(L * dx * rho* PhaseSL [i]))
            #print('Efonte',Efonte)
            #print('soustraciton',Efonte / (rho * dx * C))
            #print('Ts avant', Ts[n, i])
            Ts[n, i] = Ts[n, i]-(Efonte / (rho * dx * C))
            #print('ts apres', Ts[n, i])
            PhaseSL[i] = PhaseSL[i] - (Efonte / (L * dx * rho))
            #print('Phase fin',PhaseSL[i])
            Ts[n,i]=min(Ts[n,i],Tf)
    return Ts-273.15,PhaseSL

L=333550. #chaleur latente fiusion de la glace
rho= 917. #masse volumique
K=2.22
C=2060.

###################################################

###################################################
dx = .01  # metre   /!\ ca ne marche pas avec tous les dx si trop grand on devient absurde
Totprofond = 1
Nx = int(Totprofond / dx)
dx2 = dx * dx
Tottime = 25000
dt = 1
Nt = int(Tottime / dt)  # nb de pas de temps
print(Nt)

###################################################

epsi = 0
Tf = .0
R = 1. / (1. + 2. * C * epsi / L)

#####################################
bordhaut=10.
Tini=-2.
bordbas=-10.
#################################
#pb avec negative temperature
def Crocusfunction(L, rho, K, C, Nx, dx2, dt, Nt, Tf, bordhaut, bordbas, Tini):
    T=np.ones((Nt,Nx))*Tini
    PhaseSL=np.ones(Nx) #1=solide 0=liquide
    PhaseSL[0]=0
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
        T,PhaseSL=fonte(T, Tf, n, L, C,PhaseSL,rho,dx,Nx)
    #print(PhaseSL)
    return T
T=Crocusfunction(L, rho, K, C, Nx, dx2, dt, Nt, Tf, bordhaut, bordbas, Tini)

print("--- %s seconds ---" % (time.time() - start_time))
dir = '/Users/angehaddj/Desktop/CD/np.save/'
name = 'crocus_dx' + str(round(dx, 2)) + '*Nx' + str(round(Nx, 2)) + '*Totx' + str(
     round(Tottime, 2)) + '_dt' + str(round(dt, 2)) + '*Nt' + str(round(Nt, 2)) + '*Tott' + str(
     round(Tottime, 2)) + '_' + str(bordhaut) + str(Tini) + str(bordbas) + '.np'
np.save(dir + name, T)
Phase=phase(T,Tf,epsi)
plt.ylabel('Profondeur (m)')
plt.xlabel('Pas de Temps (s)')
xtent = [dt*0 , Nt,  dx*0, Nx]
#print(extent)
norm = mcolors.TwoSlopeNorm(vmin=T.min(), vmax = T.max(), vcenter=0) #pour fixer le 0 au blanc
im=plt.imshow(np.transpose(T),cmap='seismic' , norm=norm,aspect='auto',interpolation='None')
plt.title('CROCUS CFL: '+str(round((dt*K)/(dx2*C),5))+'\n Th: '+ str(bordhaut)+ ' Tb: '+str( bordbas)+ ' Ti: '+str(Tini)+' dt: '+str(round(dt,5))+ " dx: "+str(round(dx,5))+"\n Execution time: "+str(round(time.time() - start_time))+"s")
cb=plt.colorbar()
cb.ax.set_ylabel('Temperature °C', rotation=270)
plt.contourf(np.ma.masked_where(np.transpose(Phase)==0,(np.transpose(Phase))),0,hatches=[ '////'], alpha=0)
plt.show()
# plt.imshow(np.transpose(Phase),cmap='viridis',aspect='auto',interpolation='None')
# plt.contourf(np.ma.masked_where(np.transpose(Phase)==0,np.transpose(Phase)),0,hatches=[ '////'], alpha=0)
# plt.colorbar()
# plt.title('/!\CROCUS CFL Phase')
# plt.show()