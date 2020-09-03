# import progressbar as pb
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np


def rebin(arr, new_shape):
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)

print("load fichier")
#ref = np.load('/Users/angehaddj/Desktop/CD/np.save/verite_dx0.01*Nx200*Totx25000_dt0.01*Nt5000000*Tott25000_10.0-2.0-10.0.np.npy')
ref = np.load('/Users/angehaddj/Desktop/CD/np.save/verite_dx0.01*Nx100*Totx25000_dt0.02*Nt1250000*Tott25000_10.0-2.0-10.0.np.npy')
print('loaded')
# T=rebin(ref,(ref.shape[0]//1,ref.shape[1]//4))
# print('rebined')
# print((ref.shape[0]//1,ref.shape[1]//4))
# dir='/Users/angehaddj/Desktop/'
# name='bined3.np'
# np.save(dir+name,T)
# print(dir+name)
plt.ylabel('Profondeur (nb de pas dx)')
plt.xlabel('Temps (nb de pas dt)')
titreying = ('Schema explicit CFL: 0.48419')
print('traçage du graph')
norm = mcolors.TwoSlopeNorm(vmin=ref.min(), vmax = ref.max(), vcenter=0) #pour fixer le
# 0 au blanc
titreyang=('\n Th: 10.0 Tb: -10.0 Ti: -2.0 dt: 0.02 dx: 0.01\n Execution time: 2207s')
titre=titreying+titreyang
plt.imshow(np.transpose (ref), cmap=plt.cm.seismic, aspect='auto', interpolation='None',norm=norm)
plt.title(titre)
cb = plt.colorbar()
cb.ax.set_ylabel('Temperature °C', rotation=270)
#plt.savefig('/Users/angehaddj/Desktop/temp/'+titreying+'dt'+str(round(dt,5))+'dx'+str(round(dx,5))+'.png', format='png',dpi=400)
plt.show()
# plt.imshow(np.transpose (ref), cmap=plt.cm.seismic, aspect='auto', interpolation='None',norm=norm)
# plt.show()
# plt.imshow(np.transpose (T), cmap=plt.cm.seismic, aspect='auto', interpolation='None',norm=norm)
# plt.show()
# L = 333550.  # chaleur latente fiusion de la glace
# rho = 917.  # masse volumique
# K = 2.22
# C = 2060
# x=np.arange(0.01,0.1,0.0001)
# t=np.arange(0.01,1000,0.01)
# cfl=np.empty((x.size,t.size))
# for i in range(len(x)):
#     for j in range(len(t)):
#         cfl[i,j]=(t[j]*K)/(x[i]*x[i]*rho)
#
# norm = mcolors.TwoSlopeNorm(vmin=0, vmax = 1,vcenter=.5)
# plt.xscale('log')
# plt.ylabel('dx')
# plt.xlabel('dt')
# plt.title('Visualisation du nombre CFL \nen fonction des parametre dx et dt\nRouge: CFL>0.5')
# plt.imshow(cfl,cmap=plt.cm.seismic, aspect='auto', interpolation='None',norm=norm,extent=[0.01,1000,0,0.1])
# cb=plt.colorbar()
# cb.ax.set_ylabel('Valeur de CFL', rotation=270)
# plt.show()