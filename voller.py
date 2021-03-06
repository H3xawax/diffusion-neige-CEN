import numpy as np
import matplotlib.pyplot as plt
import progressbar as pb
import matplotlib.colors as mcolors
import time
start_time = time.time()

def Sj1 (Tj,TJ1,Tf, epsi, L, Qj1, C, R):
    if Tj <= (Tf-epsi):
        if TJ1 <= (Tf-epsi):
            #print('c\'est bon')
            return 0
        if TJ1 > (Tf+epsi): return -L
        else: return R*(-C*(epsi+Tj)-Qj1)
    if Tj >= (Tf+epsi):
        if TJ1 >= (Tf+epsi): return 0
        if TJ1 < (Tf-epsi): return L
        else: return R*(-C*(-epsi+Tj)-Qj1)
    else:
        if TJ1 > (Tf+epsi): return L*(Tj+epsi)/(2*epsi) -L
        if TJ1 <= (Tf-epsi): return L*(Tj+epsi)/(2*epsi)
        else : return -R*Qj1

def H(T,epsi, L, Tf,C):
    if T<=(Tf-epsi): return C*T
    if T>(Tf+epsi): return C*T+L
    else: return C*T+L*(T+epsi)/(2*epsi)

def dH(T,epsi, L, Tf):
    if T<=(Tf-epsi) : return 0
    if T>(Tf+epsi) : return L
    else : return L*(T+epsi)/(2*epsi)

def conv(Tj1,Tj1k, epsi, L, Tf, C,convergence):
        for i in range (len(Tj1)):
            #print('convi: ',i)
            #print('Tj1conv',Tj1[i])
            #print('Tj1kconv', Tj1k[i])
            #print('conv :', (H(Tj1[i], epsi, L, Tf, C) - H(Tj1k[i], epsi, L, Tf, C)) / C)
            convmax=0
            convcandidate=(abs(H(Tj1[i], epsi, L, Tf, C) - H(Tj1k[i], epsi, L, Tf, C)) / C)
            if convcandidate>convmax:
                convmax=convcandidate
            if (abs(H(Tj1[i], epsi, L, Tf, C) - H(Tj1k[i], epsi, L, Tf, C)) / C) > convergence :
                return False
            #else :
              #  print('Non convergence')
           # print(convmax)
        return True


def init(Tj,rho,dx2,C,dt,K,Nx):
    #Tj1[0]=Tj[0]
    Tj1=np.ones(len(Tj))
    Tj1[0]=Tj[0]
    Tj1[-1] = Tj[-1]
    #Tj1[Nx-1] = Tj[Nx-1] #pour avoir la derneiere valeur
    for i in range(1,Nx-1):
        Tj1[i]=Tj[i]+(Tj[i-1]-2*Tj[i]+Tj[i+1])*(dt*K)/(rho*dx2*C)
        #print(i,'tada:',(dt*K)/(rho*dx2*C))
    #print('Init terminated')
    #print('Tj', Tj)
    #print('Tj1',Tj1)
    return Tj,Tj1

def gaussseidel (Nx,dt,K,dx2,rho,C,epsi,L,Tf,convergence,Tj, Tj1):
    k=0
    Tj1k = np.empty_like(Tj1)
    Tj1k[:] = Tj1[:]
    lambd = 2*dt*K/(dx2*rho)
    Tjm=np.empty((Nx))
    #for f in range(50):
    #print('avant while Tj :', Tj)
    #print('avant while Tj1 :',Tj1)
    #print('avant while Tj1k :', Tj1k)
    while(True):
        for i in range(1, Nx -1):

            #print('if')
            if (0<dH(Tj1[i],epsi,L,Tf)and dH(Tj1[i],epsi,L,Tf)<L ):
                Tj1k[i]=(2*epsi*dH(Tj1[i],epsi,L,Tf))/L -epsi
                #print('if:',Tj1[i],':',Tj1k[i])



            else:

                #print('k:', k)
                #print('i:',i)
                #print('avant-Tj :', Tj[i])
                #print('avant-Tj1 :', Tj1[i])
                #print('avant-Tj1k :', Tj1k[i])
                Qj1=(Tj[i-1]-2*Tj[i]+Tj[i+1])*(dt*K)/(rho*dx2)
                #print(Tj[i],Tj[i-1],Tj[i+1])
                #print((1+lambd/C))

                Tj1k[i]= ( Tj[i] + (dt*K)*(Tj1[i-1]+Tj1[i+1])/(C*dx2*rho) + (1/C)*(Sj1(Tj[i],Tj1[i],Tf, epsi, L,Qj1 , C, R)))/(1+lambd/C)
                #if i==1:
                    #print('Sj1',Sj1(Tj[i],Tj1[i],Tf, epsi, L,Qj1 , C, R))
                    #print('Tj',Tj[i])
                    #print('Tj1', Tj1[i])
                    #print('Tj1k',Tj1k[i])
                #Tjm[i]=Tj1k[i] - Tj1[i]
            #print('Tj1k-avant: ',Tj1k[i])
            #print('Tj1-avant: ',Tj1[i])
                #print(Tj1k[i] - Tj1[i])
                #if (rho*dx2/C)*(Sj1(Tj[i],Tj1[i],Tf, epsi, L,Qj1 , C, R))!=0:print((rho*dx2/C)*(Sj1(Tj[i],Tj1[i],Tf, epsi, L,Qj1 , C, R)))
                #print(Tj[i],Tj1[i])
                #print('apres else Tj1k :', Tj1k)
            #Tj1k[i]=( Tj[i] + (dt*K)/(C*dx2*rho)*(Tj[i-1]+Tj[i+1]) + (1/C)*(dH(Tj[i],epsi,L,Tf)-dH(Tj1[i],epsi,L,Tf))) /(1+lambd/C)
                        #Tj1k[i]= Tj[i] + (1/C)*(K/(rho*dx2))*(Tj[i-1]-2*Tj[i]+Tj[i+1]) + (1/C)*(dH(Tj[i],epsi,L,Tf)-dH(Tj1[i],epsi,L,Tf))
            #print(dH(Tj[i],epsi,L,Tf)*(rho/dx2))


                #print('Tj1k-apres: ', Tj1k[i])
                #print('Tj1-apres: ', Tj1[i])
            #print('else')
            #if (conv(Tj1, Tj1k, epsi, L, Tf, C, convergence)):
            #   break
            #     return Tj1k, k
            # if (Sj1(Tj[i],Tj1[i],Tf, epsi, L,Qj1 , C, R)!=0):
            #     print(Tj1k[i])
            #     Tj1k[i]=(2*epsi*dH(Tj1[i],epsi,L,Tf))/L -epsi
            #     print('if:',Tj1[i],':',Tj1k[i])
        #print('fin de boucle Tj : ',Tj)
        #print('fin de boucle Tj1 : ',Tj1)
        #print('fin de boucle Tj1k : ',Tj1k)
        #print('k:',k)
        #print('while:',(H(Tj1[i], epsi, L, Tf, C) - H(Tj1k[i], epsi, L, Tf, C)*(rho*dx2/C)))
        k = k + 1
        #print(np.mean(Tjm))
        #print('k:', k)
        #for i in range (len(Tj1)):
            # print('boucle: ',i)
            # print('Tj1k',Tj1k[i])
            # print('Tj1',Tj1[i])
            # print('diff :', (H(Tj1[i], epsi, L, Tf, C) - H(Tj1k[i], epsi, L, Tf, C)) / C)
        if (conv(Tj1,Tj1k, epsi, L, Tf, C,convergence)): break
        #if k>=100:break #pour eviter les non convergence
        Tj1[:]=Tj1k[:]
            #print('Sj+1',Sj1(Tj[i], Tj1[i], Tf, epsi, L, Qj1, C, R))


    return Tj1k,k
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
Tottime=500000
dt=10000# seconde
Nt= int(Tottime/dt)#nb de pas de temps
Totprofond=1
dx=.01#np.sqrt(dt*K/(C*.43)) #  metre   /!\ ca ne marche pas avec tous les dx si trop grand on devient absurde
Nx=int(Totprofond/dx)
dx2=dx*dx
print(dt)
print(dx2)
print(Nt)
print(Nx)


###################################################
# dx=.05 #  metre   /!\ ca ne marche pas avec tous les dx si trop grand on devient absurde
# Totprofond=1
# Nx=int(Totprofond/dx)
# dx2=dx*dx
# Tottime=50000
# dt=dx2*C*.12/K
# Nt= int(Tottime/dt)#nb de pas de temps

###################################################

epsi=.5
Tf=.0
convergence=.001
R=1./(1.+2.*C*epsi/L)
lambd = 2*dt*K/(dx2*rho)
#print(dx)
#print((1+lambd/C))
#print(L/(C*(1+lambd/C)))

print((L/C)/(2+(2*lambd/C)))

#####################################
bordhaut=2.
Tini=2.
bordbas=-5.
#################################


T=np.ones((Nt,Nx))*Tini
#T=np.ones((Nt,Nx))*np.linspace(bordhaut-5,bordbas,Nx)

T[:,0]=bordhaut #voir linspace
T[:,Nx-1]=bordbas
kmax=0
print('############### CFL: ',(dt*K)/(dx2*C),"###############")

for t in pb.progressbar(range(1,Nt)):
    #print('######t: ', t)
    Tj,Tj1=init(T[t - 1, :], rho, dx2, C, dt, K, Nx)
    T[t,:],kcandidate=gaussseidel(Nx,dt,K,dx2,rho,C,epsi,L,Tf,convergence,Tj,Tj1)
    if kcandidate>kmax:
        kmax=kcandidate
        print('kmax: ',kmax)
    #if t>=4350:print(T[t,:])
    #junk,T[t,:]=init(T[t-1,:],rho,dx2,C, dt, K, Nx)
print(':', T[Nt-1,:])


print("--- %s seconds ---" % (time.time() - start_time))
Phase=phase(T,Tf,epsi)

plt.ylabel('Profondeur (m)')
plt.xlabel('Temps (s)')
extent = [dt*0 , Nt,  dx*0, Nx]
#print(extent)
norm = mcolors.TwoSlopeNorm(vmin=T.min(), vmax = T.max(), vcenter=0) #pour fixer le 0 au blanc
im=plt.imshow(np.transpose(T),cmap=plt.cm.seismic, norm=norm ,aspect='auto',interpolation='None')
plt.title('VOLLER CFL: '+str(round((dt*K)/(dx2*C),5))+'\n Th: '+ str(bordhaut)+ ' Tb: '+str( bordbas)+ ' Ti: '+str(Tini)+' dt: '+str(round(dt,5))+ " dx: "+str(round(dx,5))+"\n Execution time: "+str(round(time.time() - start_time))+"s")
cb=plt.colorbar()
cb.ax.set_ylabel('Temperature °C', rotation=270)
plt.contourf(np.transpose(Phase),5,hatches=[ '','////','////','////','////'], alpha=0,aspect='auto',interpolation='None')
plt.show()
