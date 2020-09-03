import time

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import progressbar as pb

start_time = time.time()


def conv(Tj1, Tj1k, epsi, L, Tf, C, convergence):
    for i in range(len(Tj1)):
        # print('convi: ',i)
        # print('Tj1conv',Tj1[i])
        # print('Tj1kconv', Tj1k[i])
        # print('conv :', (H(Tj1[i], epsi, L, Tf, C) - H(Tj1k[i], epsi, L, Tf, C)) / C)
        convmax = 0
        convcandidate = (abs(Phi_1(Tj1[i], epsi, L, Tf, C) - Phi_1(Tj1k[i], epsi, L, Tf, C)) / C)
        if (convcandidate > convmax):
            convmax = convcandidate
        if abs(Phi_1(Tj1[i], epsi, L, Tf, C) - Phi_1(Tj1k[i], epsi, L, Tf, C)) / C > convergence:
            return False
        # else :
        #  print('Non convergence')
    # print(convmax)
    return True

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
dx = .05  # metre   /!\ ca ne marche pas avec tous les dx si trop grand on devient absurde
Totprofond = 1
Nx = int(Totprofond / dx)
dx2 = dx * dx
Tottime = 25000
dt = 1
Nt = int(Tottime / dt)  # nb de pas de temps
print(Nt)

###################################################
##@@@@>
###############Y A MOYEN IL MANQUE UN RHO DANS LE CFL dx2*C*.45*rho/K
############http://www-udc.ig.utexas.edu/external/becker/teaching/557/problem_sets/problem_set_fd_implicit.pdf

###################################################
epsi=.1
Tf=.0
lambd = 2*dt*K/(dx2*rho)
#print(dx)
#print((1+lambd/C))
#print(L/(C*(1+lambd/C)))

#print((L/C)/(2+(2*lambd/C)))

#####################################
bordhaut=10.
Tini=-2.
bordbas=-10.
#################################

def voller2function(L,rho,K,C,Nx,dx2,dt,Nt,epsi,Tf,bordhaut,bordbas,Tini,lambd):
    T=np.ones((Nt,Nx))*Tini
    #T=np.ones((Nt,Nx))*np.linspace(bordhaut-5,bordbas,Nx)

    T[:,0]=bordhaut #voir linspace
    T[:,Nx-1]=bordbas
    print('############### CFL: ',(dt*K)/(dx2*C),"###############")

    for t in pb.progressbar(range(Nt-1)):
        while (True):
            T1=T
            k = 0
            for i in range(1,Nx-1):
                T1[t+1,i]=( Phi_1(T[t,i],Tf, epsi,L,C) + (dt*K)*(T[t,i-1]+T[t,i+1])/(dx2*rho) + alpha(T[t,i], Tf, epsi,L))/beta(T[t,i],Tf, epsi,L,C,lambd)
                k = k + 1
                print(k)
            if conv(T, T1, t, .01, Nx): break
            T = T1
    return T

def conv (T,T1,t,convergence,Nx):
    for i in range(1,Nx-1) :
        if abs(Phi_1(T[t,i], Tf, epsi,L,C)- Phi_1(T1[t,i], Tf, epsi,L,C))/C > convergence:
            return False
    return True


T=voller2function(L,rho,K,C,Nx,dx2,dt,Nt,epsi,Tf,bordhaut,bordbas,Tini,lambd)
print("--- %s seconds ---" % (time.time() - start_time))

Phase=phase(T,Tf,epsi)
print(np.shape(T))
print(np.shape(Phase))
# dir = '/Users/angehaddj/Desktop/CD/np.save/'
# name = 'verite_voller_dx' + str(round(dx, 2)) + '*Nx' + str(round(Nx, 2)) + '*Totx' + str(
#     round(Tottime, 2)) + '_dt' + str(round(dt, 2)) + '*Nt' + str(round(Nt, 2)) + '*Tott' + str(
#     round(Tottime, 2)) + '_' + str(bordhaut) + str(Tini) + str(bordbas) + '.np'
# np.save(dir + name, T)
plt.ylabel('Profondeur (m)')
plt.xlabel('Pas de temps (s)')
extent = [dt*0 , Nt,  dx*0, Nx]
#print(extent)
norm = mcolors.TwoSlopeNorm(vmin=-10, vmax = 10, vcenter=0) #pour fixer le 0 au blanc
im=plt.imshow(np.transpose(T),cmap=plt.cm.seismic, norm=norm ,aspect='auto',interpolation='None',extent=extent)
plt.title('/!\VOLLER2 CFL: '+str(round((dt*K)/(dx2*C),5))+'\n Th: '+ str(bordhaut)+ ' Tb: '+str( bordbas)+ ' Ti: '+str(Tini)+' dt: '+str(round(dt,5))+ " dx: "+str(round(dx,5))+"\n Hatching is liquid phase\n Execution time: "+str(round(time.time() - start_time))+"s")
cb=plt.colorbar()
cb.ax.set_ylabel('Temperature °C', rotation=270)
#plt.contourf(np.ma.masked_where(np.transpose(Phase)==0,(np.transpose(Phase))),0,hatches=[ '////'], alpha=0)
plt.show()
# plt.imshow(np.transpose(Phase),cmap='Greys',aspect='auto',interpolation='None')
# plt.colorbar()
# plt.title('/!\VOLLER2 CFL Phase')
# plt.show()