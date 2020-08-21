import numpy as np
import matplotlib.pyplot as plt
# import progressbar as pb
import matplotlib.colors as mcolors
import time
from fonction_comparaison import phase, explicitfunction, voller2function

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
dt = 625 # multiple 1 2 4 5 8 10 20 25 40 50 100 125 200 250 500 625 1000
Nt = int(Tottime / dt)  # nb de pas de temps
print(Nt)

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
lambd = 2 * dt * K / (dx2 * rho)

# T est transposé dans ce qui suit !!!!!!!!!!!!


######################################
print('ENTREZ un numero\n 1 explicit\n 2 voller2\n 99 importer un fichier')
select = None
while not select:
    try:
        select = int(input('>'))
    except ValueError:
        print('Invalid Number')
if select == 1:
    T = explicitfunction(L, rho, K, C, Nx, dx2, dt, Nt, epsi, Tf, bordhaut, bordbas, Tini)
if select == 2:
    T = voller2function(L, rho, K, C, Nx, dx2, dt, Nt, epsi, Tf, bordhaut, bordbas, Tini, lambd)
if select == 99:
    T = np.load(
        '/Users/angehaddj/Desktop/CD/np.save/verite_voller_dx0.01*Nx100*Totx25000_dt1000*Nt25*Tott25000_10.0-2.0-10.0.np.npy')
print("--- %s seconds ---" % (time.time() - start_time))

Phase = phase(T, Tf, epsi)
if select == 1:
    titre = 'EXPLICIT CFL: '
if select == 2:
    titre = 'VOLLER2 CFL: '
#
plt.ylabel('Profondeur (nb de pas)')
plt.xlabel('Temps (nb de pas)')
# xtent = [dt*0 , Nt,  dx*0, Nx]
# #print(extent)
norm = mcolors.TwoSlopeNorm(vcenter=0)  # pour fixer le 0 au blanc
# im=plt.imshow(np.transpose(T),cmap=plt.cm.seismic, norm=norm ,aspect='auto',interpolation='None')
# plt.title(titre+str(round((dt*K)/(dx2*C),5))+'\n Th: '+ str(bordhaut)+ ' Tb: '+str( bordbas)+ ' Ti: '+str(Tini)+' dt: '+str(round(dt,5))+ " dx: "+str(round(dx,5))+"\n Execution time: "+str(round(time.time() - start_time))+"s")
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

save = np.load('/Users/angehaddj/Desktop/CD/np.save/verite_dx0.01*Nx100*Totx25000_dt0.04*Nt625000*Tott25000_10.0-2.0-10.0.np.npy')
print(np.shape(T)[0])
#####################################################################################################
#          reformat pour la soustraction comparaison
print(((np.shape(save)[0] / np.shape(T)[0]), (np.shape(save)[1] / np.shape(T)[1])))
T = np.kron(T, np.ones((int(np.shape(save)[0] / np.shape(T)[0]), int(np.shape(save)[1] / np.shape(T)[1]))))
#####################################################################################################


plt.imshow(np.transpose(save - T), cmap=plt.cm.seismic, aspect='auto', interpolation='None',norm=norm)
plt.title('save-calcul')
cb = plt.colorbar()
cb.ax.set_ylabel('Temperature °C', rotation=270)
plt.show()
