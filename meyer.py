import numpy as np
import matplotlib.pyplot as plt
import progressbar as pb
import matplotlib.colors as mcolors
import time
start_time = time.time()

def H(u,Tf,epsi,ci,cw,L):
    if u<=Tf-epsi:return ci*u
    if u>Tf+epsi:return H(Tf+epsi,Tf,epsi,ci,cw,L)+cw*(u-Tf-epsi)
    else:return H(Tf-epsi,Tf,epsi,ci,cw,L)+L*(u-Tf-epsi)/(2*epsi)

def Vi(i,t,T,K,rho,dx2,dt,force,Tf,epsi,ci,cw,L):
    if force==1:
        return  (H(T[t-1,i],Tf,epsi,ci,cw,L)/dt -K*(T[t,i-1]+T[t-1,i+1])/(rho*dx2))/(-2/(rho*dx2)+ci/dt)
    if force==3:
        return  (cw*(Tf+epsi)/dt-H(Tf+epsi,Tf,epsi,ci,cw,L)/dt+H(T[t-1,i],Tf,epsi,ci,cw,L)/dt -K*(T[t,i-1]+T[t-1,i+1])/(rho*dx2))/(-2/(rho*dx2)+cw/dt)
    if force==2:
        return  (L*(Tf+epsi)/(2*epsi*dt)-H(Tf-epsi,Tf,epsi,ci,cw,L)/dt+H(T[t-1,i],Tf,epsi,ci,cw,L)/dt -K*(T[t,i-1]+T[t-1,i+1])/(rho*dx2))/(-2/(rho*dx2)+L/(dt*2*epsi))
    else: print("### FORCE ERROR ###")

def verifvi(vi,Tf,epsi):
    if vi<=Tf-epsi: return 1
    if vi>=Tf+epsi: return 3
    else: return 2
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
ci=2060
cw=4185
###################################################
dx=.025 #  metre   /!\ ca ne marche pas avec tous les dx si trop grand on devient absurde
Totprofond=1
Nx=int(Totprofond/dx)
dx2=dx*dx
Tottime=5000
dt=dx2*ci*.45/K
Nt= int(Tottime/dt)#nb de pas de temps
###################################################
epsi=.0001
Tf=.0
lambd = 2*dt*K/(dx2*rho)
#####################################
bordhaut=20.
Tini=-2.
bordbas=20.
#################################
convergence=0.0001
T=np.ones((Nt,Nx))*Tini

T[:,0]=bordhaut #voir linspace
T[:,Nx-1]=bordbas
print('############### CFL: ',(dt*K)/(dx2*ci),"###############")
A=np.zeros((Nx,Nx))
for i in range(1,Nx-1):
    A[i,i-1]=K/(rho*dx2)
    A[i,i]=-2*K/(rho*dx2)
    A[i, i + 1] = K / (rho * dx2)

SSORparam=.1
print('############### SSORparam: ',SSORparam,"###############")

for t in pb.progressbar(range(1,Nt-1)):
    for i in range(1,Nx-1):
        force1=1#verifvi(T[t-1,i],Tf,epsi)
        vi=Vi(i, t, T, K, rho, dx2, dt, force1, Tf, epsi, ci, cw, L)
        force2 = verifvi(vi, Tf, epsi)
        while(force1!=force2):
            vi = Vi(i, t, T, K, rho, dx2, dt, force2, Tf, epsi, ci, cw, L)
            force3=force2
            force2=verifvi(vi,Tf,epsi)
            force1=force3
        temp1=T[t-1,i]+SSORparam*(vi-T[t-1,i])
        temp2=2000#grande valeur pour l'initialisation du while
        while(True):
            #print("while")
            force1 = verifvi(T[t - 1, i], Tf, epsi)
            vi = Vi(i, t, T, K, rho, dx2, dt, force1, Tf, epsi, ci, cw, L)
            force2 = verifvi(vi, Tf, epsi)
            while (force1 != force2):
                vi = Vi(i, t, T, K, rho, dx2, dt, force2, Tf, epsi, ci, cw, L)
                force3 = force2
                force2 = verifvi(vi, Tf, epsi)
                force1 = force3
            temp2 = T[t - 1, i] + SSORparam * (vi - T[t - 1, i])
            if (abs(temp1-temp2)<convergence):break
            temp1=temp2



print("--- %s seconds ---" % (time.time() - start_time))
Phase=phase(T,Tf,epsi)

plt.ylabel('Profondeur (m)')
plt.xlabel('Temps (s)')
extent = [dt*0 , Nt,  dx*0, Nx]
#print(extent)
norm = mcolors.TwoSlopeNorm(vmin=T.min(), vmax = T.max(), vcenter=0) #pour fixer le 0 au blanc
im=plt.imshow(np.transpose(T),cmap=plt.cm.seismic, norm=norm ,aspect='auto',interpolation='None')
plt.title('WHITE CFL: '+str(round((dt*K)/(dx2*ci),5))+'\n Th: '+ str(bordhaut)+ ' Tb: '+str( bordbas)+ ' Ti: '+str(Tini)+' dt: '+str(round(dt,5))+ " dx: "+str(round(dx,5))+"\n Execution time: "+str(round(time.time() - start_time))+"s")
cb=plt.colorbar()
cb.ax.set_ylabel('Temperature °C', rotation=270)
plt.contourf(np.transpose(Phase),1,hatches=['', '////'], alpha=0,aspect='auto',interpolation='None')
plt.show()
