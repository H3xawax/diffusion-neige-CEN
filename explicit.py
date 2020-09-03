import time

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import progressbar as pb

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
K=2.22#conductivité
C=2060

###################################################
dx = .1  # metre   /!\ ca ne marche pas avec tous les dx si trop grand on devient absurde
Totprofond = 1
Nx = int(Totprofond / dx)
dx2 = dx * dx
Tottime = 25000
dt = 2
Nt = int(Tottime / dt)  # nb de pas de temps
print(Nt)


###################################################

epsi=.0
Tf=.0
#convergence=.001
#R=1./(1.+2.*C*epsi/L)

###################################################
bordhaut=10.
Tini=-2.
bordbas=-10.
###################################################

def explicitfunction(L,rho,K,C,Nx,dx2,dt,Nt,epsi,Tf,bordhaut,bordbas,Tini):
    T=np.ones((Nt,Nx))*Tini
    H=np.empty_like(T)
    #T=np.ones((Nt,Nx))*np.linspace(bordhaut-5,bordbas,Nx)

    T[:,0]=bordhaut #voir linspace
    T[:,Nx-1]=bordbas
    print('############### CFL: ',(dt*K)/(dx2*rho),"###############")
    print('Transformation de la matrice en Entalpie')
    for t in pb.progressbar(range(Nt)):
        for x in range(Nx):
            H[t,x] = Hconv(T[t,x],epsi, L, Tf,C)
    kmax=0
    print()
    print('############### CFL: ',(dt*K)/(dx2*rho),"###############")
    print('Résolution de la diffusion en Entalpie')

    for t in pb.progressbar(range(1,Nt)):
        #print('######t: ', t)
        H[t,:]=explicit(H[t-1,:],rho,dx2,C,dt,K,Nx)

    print('Transformation de la matrice en Temperature')
    for t in pb.progressbar(range(Nt)):
        for x in range(Nx):
            T[t,x] = Tconv(H[t,x],epsi, L,C)
    #print(':', T[Nt-1,:])
    return T
T=explicitfunction(L,rho,K,C,Nx,dx2,dt,Nt,epsi,Tf,bordhaut,bordbas,Tini)
print("--- %s seconds ---" % (time.time() - start_time))
Phase=phase(T,Tf,epsi)

dir='/Users/angehaddj/Desktop/CD/np.save/'
name='verite_dx'+str(round(dx,2))+'*Nx'+str(round(Nx,2))+'*Totx'+str(round(Tottime,2))+'_dt'+str(round(dt,2))+'*Nt'+str(round(Nt,2))+'*Tott'+str(round(Tottime,2))+'_'+str(bordhaut)+str(Tini)+str(bordbas)+'.np'
np.save(dir+name,T)
print(dir+name)
#
# plt.ylabel('Nb de pas de profondeur')
# plt.xlabel('Nb de pas de temps')
# #xtent = [dt*0 , Nt,  dx*0, Nx]
# #print(extent)
# norm = mcolors.TwoSlopeNorm(vmin=T.min(), vmax = T.max(), vcenter=10) #pour fixer le 0 au blanc
# im=plt.imshow(np.transpose(T),cmap=plt.cm.seismic, norm=norm ,aspect='auto',interpolation='None')
# plt.title('EXPLICIT CFL: '+str(round((dt*K)/(dx2*C),5))+'\n Th: '+ str(bordhaut)+ ' Tb: '+str( bordbas)+ ' Ti: '+str(Tini)+' dt: '+str(round(dt,5))+ " dx: "+str(round(dx,5))+"\n Execution time: "+str(round(time.time() - start_time))+"s")
# cb=plt.colorbar()
# cb.ax.set_ylabel('Temperature °C', rotation=270)
# #plt.contourf(np.ma.masked_where(np.transpose(Phase)==0,np.transpose(Phase)),0,hatches=[ '////'], alpha=0)
# plt.show()
# # plt.imshow(np.transpose(Phase),cmap='viridis',aspect='auto',interpolation='None')
# # plt.contourf(np.ma.masked_where(np.transpose(Phase)==0,np.transpose(Phase)),0,hatches=[ '////'], alpha=0)
# # plt.colorbar()
# # plt.title('/!\explicit CFL Phase')
# # plt.show()
plt.ylabel('Profondeur (Nb de pas)')
plt.xlabel('Pas de temps (nb de pas)')
extent = [dt*0 , Nt,  dx*0, Nx]
#print(extent)
norm = mcolors.TwoSlopeNorm(vmin=-10, vmax = 10, vcenter=0) #pour fixer le 0 au blanc
im=plt.imshow(np.transpose(T),cmap=plt.cm.seismic, norm=norm ,aspect='auto',interpolation='None',extent=extent)
plt.title('/!\explicit CFL: '+str(round((dt*K)/(dx2*rho),5))+'\n Th: '+ str(bordhaut)+ ' Tb: '+str( bordbas)+ ' Ti: '+str(Tini)+' dt: '+str(round(dt,5))+ " dx: "+str(round(dx,5))+"\n Execution time: "+str(round(time.time() - start_time))+"s")
cb=plt.colorbar()
cb.ax.set_ylabel('Temperature °C', rotation=270)
#plt.contourf(np.ma.masked_where(np.transpose(Phase)==0,np.transpose(Phase)),0,hatches=[ '////'], alpha=0)
plt.show()
# plt.imshow(np.transpose(Phase),cmap='Greys',aspect='auto',interpolation='None')
# plt.colorbar()
# plt.title('/!\VOLLER2 CFL Phase')
# plt.show()