import time

# import progressbar as pb
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np

from fonction_comparaison import phase, explicitfunction, voller2function, Crocusfunction

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
dt = 1# multiple 1 2 4 5 8 10 20 25 40 50 100 125 200 250 500 625 1000
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
print('ENTREZ un numero\n 1 Explicit\n 2 Voller2\n 3 Crocus\n 99 importer un fichier')
select = None
T=None
titreying=None
while not select:
    try:
        select = int(input('>'))
    except ValueError:
        print('Invalid Number')
if select == 1:
    titreying='Comparaison explicit CFL: '
    T = explicitfunction(L, rho, K, C, Nx, dx2, dt, Nt, epsi, Tf, bordhaut, bordbas, Tini)
if select == 2:
    titreying='Comparaison voller CFL: '
    T = voller2function(L, rho, K, C, Nx, dx2, dt, Nt, epsi, Tf, bordhaut, bordbas, Tini, lambd,convergence)
if select == 3:
    titreying='Comparaison Crocus CFL: '
    T = Crocusfunction(L, rho, K, C, Nx,dx, dx2, dt, Nt, Tf, bordhaut, bordbas, Tini)
if select == 99:
    titreying="Voir le fichier loadé"
    T = np.load(
        '/Users/angehaddj/Desktop/CD/np.save/verite_dx0.05*Nx20*Totx25000_dt1*Nt25000*Tott25000_10.0-2.0-10.0.np.npy')
        #'/Users/angehaddj/Desktop/CD/np.save/verite_dx0.05*Nx20*Totx25000_dt1*Nt25000*Tott25000_10.0-2.0-10.0.np.npy')
print("--- %s seconds ---" % (time.time() - start_time))

Phase = phase(T, Tf, epsi)

#
plt.ylabel('Profondeur (nb de pas)')
plt.xlabel('Temps (nb de pas)')
# extent = [dt*0 , Nt,  dx*0, Nx]
#print(extent)
norm = mcolors.TwoSlopeNorm(vcenter=0)  # pour fixer le 0 au blanc
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
print("load fichier")
ref = np.load('/Users/angehaddj/Desktop/CD/np.save/verite_dx0.01*Nx100*Totx25000_dt0.04*Nt625000*Tott25000_10.0-2.0-10.0.np.npy')
print(np.shape(T)[0])

print("reformat")
#####################################################################################################
#          reformat pour la soustraction comparaison
print(((np.shape(ref)[0] / np.shape(T)[0]), (np.shape(ref)[1] / np.shape(T)[1])))
T = np.kron(T, np.ones((int(np.shape(ref)[0] / np.shape(T)[0]), int(np.shape(ref)[1] / np.shape(T)[1]))))
#####################################################################################################
print(np.std(abs(T-ref)))

print('traçage du graph')
norm = mcolors.TwoSlopeNorm(vmin=T.min(), vmax = T.max(), vcenter=0) #pour fixer le
# 0 au blanc
titreyang=str(round((dt*K)/(dx2*C),5))+'\n Th: '+ str(bordhaut)+ ' Tb: '+str( bordbas)+ ' Ti: '+str(Tini)+' dt: '+str(round(dt,5))+ " dx: "+str(round(dx,5))+"\n Execution time: "+str(round(time.time() - start_time))+"s\n Note: (T-Tref)>0 <=> T>Tref"
titre=titreying+titreyang
plt.imshow(np.transpose (T-ref), cmap=plt.cm.seismic, aspect='auto', interpolation='None',norm=norm)
plt.title(titre)
cb = plt.colorbar()
cb.ax.set_ylabel('Temperature °C', rotation=270)
plt.savefig('/Users/angehaddj/Desktop/temp/'+titreying+'dt'+str(round(dt,5))+'dx'+str(round(dx,5))+'.png', format='png',dpi=1200)
plt.show()
plt.imshow(np.transpose (ref), cmap=plt.cm.seismic, aspect='auto', interpolation='None',norm=norm)
plt.show()
plt.imshow(np.transpose (T), cmap=plt.cm.seismic, aspect='auto', interpolation='None',norm=norm)
plt.show()