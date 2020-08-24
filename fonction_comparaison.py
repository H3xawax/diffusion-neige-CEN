import numpy as np
import matplotlib.pyplot as plt
import progressbar as pb
import matplotlib.colors as mcolors
import time

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

#############################################################################
#8888888888                   888 d8b          d8b 888
#888                          888 Y8P          Y8P 888
#888                          888                  888
#8888888    888  888 88888b.  888 888  .d8888b 888 888888
#888        `Y8bd8P' 888 "88b 888 888 d88P"    888 888
#888          X88K   888  888 888 888 888      888 888
#888        .d8""8b. 888 d88P 888 888 Y88b.    888 Y88b.
#8888888888 888  888 88888P"  888 888  "Y8888P 888  "Y888
#                    888
#                    888
#                    888
#############################################################################
def explicitfunction(L,rho,K,C,Nx,dx2,dt,Nt,epsi,Tf,bordhaut,bordbas,Tini):
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
        H[t,:]=explicit(H[t-1,:],rho,dx2,C,dt,K,Nx,epsi,L)

    print('Transformation de la matrice en Temperature')
    for t in pb.progressbar(range(Nt)):
        for x in range(Nx):
            T[t,x] = Tconv(H[t,x],epsi, L,C)
    #print(':', T[Nt-1,:])
    return T

def Hconv(T,epsi, L, Tf,C):
    if T<=(Tf-epsi): return C*T
    if T>(Tf+epsi): return C*T+L
    else: return C*T+L*(T+epsi)/(2*epsi)

def Tconv(H,epsi,L,C):
    if H<=(-C*epsi) : return H/C
    if H>(C*epsi+L) : return (H-L)/C
    else : return epsi*(H-L/2)/(C*epsi+L/2)


def explicit(Hj,rho,dx2,C,dt,K,Nx,epsi,L):
    Hj1=np.ones(len(Hj))
    Hj1[0]=Hj[0]
    Hj1[-1] = Hj[-1]
    for i in range(1,Nx-1):
        Hj1[i]=Hj[i]+(Tconv(Hj[i-1],epsi,L,C)-2*Tconv(Hj[i],epsi,L,C)+Tconv(Hj[i+1],epsi,L,C))*(dt*K)/(rho*dx2)
        #print(i,'tada:',(dt*K)/(rho*dx2*C))
    return Hj1
#############################################################################
#888     888          888 888
#888     888          888 888
#888     888          888 888
#Y88b   d88P  .d88b.  888 888  .d88b.  888d888
# Y88b d88P  d88""88b 888 888 d8P  Y8b 888P"
#  Y88o88P   888  888 888 888 88888888 888
#   Y888P    Y88..88P 888 888 Y8b.     888
#    Y8P      "Y88P"  888 888  "Y8888  888
#############################################################################

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

def voller2function(L,rho,K,C,Nx,dx2,dt,Nt,epsi,Tf,bordhaut,bordbas,Tini,lambd):
    T=np.ones((Nt,Nx))*Tini
    #T=np.ones((Nt,Nx))*np.linspace(bordhaut-5,bordbas,Nx)

    T[:,0]=bordhaut #voir linspace
    T[:,Nx-1]=bordbas
    print('############### CFL: ',(dt*K)/(dx2*C),"###############")

    for t in pb.progressbar(range(Nt-1)):
        for i in range(1,Nx-1):
            T[t+1,i]=( Phi_1(T[t,i],Tf, epsi,L,C) + (dt*K)*(T[t,i-1]+T[t,i+1])/(dx2*rho) + alpha(T[t,i], Tf, epsi,L))/beta(T[t,i],Tf, epsi,L,C,lambd)
    return T


##############################################
###################CROCUS#####################
##############################################

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

def Crocusfunction(L, rho, K, C, Nx, dx,dx2, dt, Nt, Tf, bordhaut, bordbas, Tini):
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