# import progressbar as pb
import time

import matplotlib.pyplot as plt
import numpy as np

from fonction_comparaison import explicitfunction

start_time = time.time()

L = 333550.  # chaleur latente fiusion de la glace
rho = 917.  # masse volumique
K = 2.22
C = 2060

###################################################
dx = .01  # metre   /!\ ca ne marche pas avec tous les dx si trop grand on devient absurde
Totprofond = 1
Nx = int(Totprofond / dx)
dx2 = dx * dx
Tottime = 25000


###################################################
# if float(625000/Nt).is_integer()==False:
#     print('Mince')
#     dt=float(625000/(int(625000/Nt)*Tottime))
#     print('Le nouveau dt=',dt)
###################################################

epsi = 1
Tf = .0
convergence = .001
R = 1. / (1. + 2. * C * epsi / L)

#####################################
bordhaut = 10.
Tini = -2.
bordbas = -10.
#################################
print("load fichier")
ref = np.load('/Users/angehaddj/Desktop/CD/np.save/verite_dx0.01*Nx100*Totx25000_dt0.04*Nt625000*Tott25000_10.0-2.0-10.0.np.npy')
def compa(L, rho, K, C, Nx, dx, dx2, dt, Nt, Tf, bordhaut, bordbas, Tini,ref):
    lambd = 2 * dt * K / (dx2 * rho)
    #T=Crocusfunction(L, rho, K, C, Nx, dx, dx2, dt, Nt, Tf, bordhaut, bordbas, Tini)
    #T=voller2function(L, rho, K, C, Nx, dx, dx2, dt, Nt, Tf, bordhaut, bordbas, Tini,lambd,0.001)
    T=explicitfunction(L, rho, K, C, Nx, dx, dx2, dt, Nt, Tf, bordhaut, bordbas, Tini)

    # T est transposé dans ce qui suit !!!!!!!!!!!!
    print("--- %s seconds ---" % (time.time() - start_time))

    #Phase = phase(T, Tf, epsi)

    #

    # extent = [dt*0 , Nt,  dx*0, Nx]
    #print(extent)
    # im=plt.imshow(np.transpose(T),cmap=plt.cm.seismic, norm=norm ,aspect='auto',interpolation='None')
    # plt.title(titreying+str(round((dt*K)/(dx2*C),5))+'\n Th: '+ str(bordhaut)+ ' Tb: '+str( bordbas)+ ' Ti: '+str(Tini)+' dt: '+str(round(dt,5))+ " dx: "+str(round(dx,5))+"\n Execution time: "+str(round(time.time() - start_time))+"s")
    # cb=plt.colorbar()
    # cb.ax.set_ylabel('Temperature °C', rotation=270)
    # plt.contourf(np.ma.masked_where(np.transpose(Phase)==0,np.transpose(Phase)),0,hatches=[ '////'], alpha=0)
    # plt.show()

    # plt.imshow(np.transpose(Phase),cmap='viridis',aspect='auto',interpolation='None')
    # plt.contourf(np.ma.masked_where(np.transpose(Phase)==0,np.transpose(Phase)),0,hatches=[ '////'], alpha=0)
    # plt.colorbar()
    # plt.title('/!\comparaison CFL Phase')
    # plt.show()

    # dir = '/Users/angehaddj/Desktop/CD/np.save/'
    # name = 'verite_voller_dx' + str(round(dx, 2)) + '*Nx' + str(round(Nx, 2)) + '*Totx' + str(
    #     round(Tottime, 2)) + '_dt' + str(round(dt, 2)) + '*Nt' + str(round(Nt, 2)) + '*Tott' + str(
    #     round(Tottime, 2)) + '_' + str(bordhaut) + str(Tini) + str(bordbas) + '.np'
    # np.save(dir + name, T)


    print("reformat")
    #####################################################################################################
    #          reformat pour la soustraction comparaison
    print(((np.shape(ref)[0] / np.shape(T)[0]), (np.shape(ref)[1] / np.shape(T)[1])))
    T = np.kron(T, np.ones((int(np.shape(ref)[0] / np.shape(T)[0]), int(np.shape(ref)[1] / np.shape(T)[1]))))
    #####################################################################################################
    return T-ref

# T=np.empty((3,625000))
dt=np.array([1,10,100,500,1000])
for i in range(len(dt)):
    #print('ok')
    #print(dt[i])
    Nt = int(Tottime / dt[i])  # nb de pas de temps
    T=(compa(L, rho, K, C, Nx, dx, dx2, dt[i], Nt, Tf, bordhaut, bordbas, Tini,ref)[-1,:])
    #print('ok')
    #print(np.shape(T))
    print("Plot")
    plt.plot(T,np.arange(100),label=str(dt[i])+'s ')
    #plt.gca().invert_yaxis()

plt.ylabel('Profondeur (nb de pas)')
plt.xlabel('Erreur par / ref  (°C)')
plt.title("Erreur du gradient de temperature\ndans la profondeur a la fin de simulation (25000s)\nen fonction des pas de temps de simulation")
# plt.axvline(x=0.,color='k',ls=':')
plt.legend()
#plt.gca().invert_yaxis()
plt.grid(alpha=.3)
plt.gca().invert_yaxis()
plt.show()