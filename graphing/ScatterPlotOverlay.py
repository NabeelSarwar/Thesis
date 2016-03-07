import matplotlib
import gc
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.markers as markers
plt.ioff()

import numpy as np
import pyfits
matchedCat = pyfits.getdata('MatchHellDeconv2.fits')

starindices = np.where(matchedCat['mu_class'] == 2)[0]
galaxyindices = np.where(matchedCat['mu_class'] == 1)[0]

print starindices
print galaxyindices
starGalHash = {}
for e in starindices:
    starGalHash[e] = 'star'

for e in galaxyindices:
    starGalHash[e] = 'galaxy'

magG = matchedCat['magG']
magR = matchedCat['magR']
magI = matchedCat['magI']
magZ = matchedCat['magZ']
magY = matchedCat['magY']

ymags = matchedCat['mag_y']
jmags = matchedCat['mag_j']
hmags = matchedCat['mag_h']
kmags = matchedCat['mag_k']

chan1mags = matchedCat['mag_36']
chan2mags = matchedCat['mag_45']
chan3mags = matchedCat['mag_58']
chan4mags = matchedCat['mag_80']

def randomPlot(magSource1, magSource2):
    goodIndices = np.array([True for i in range(len(magSource1))])
    numbers = np.random.permutation(np.where(goodIndices == True)[0])
    print len(numbers)
    for e in numbers:
        if starGalHash[e] == 'star':
            plt.plot(magSource1[e], magSource2[e], markersize=1, marker='.', c='red')
        else:
            plt.plot(magSource1[e], magSource2[e], markersize=1, marker='.', c='blue')
    # speed up performance
    gc.collect()
    print 'done'

plt.figure()
ax1= plt.axes()
gr = magG-magR
ri = magR-magI
ax1.set_xlabel('g-r')
ax1.set_ylabel('r-i')
randomPlot(gr, ri)
plt.savefig('data/overlay/grvsri.png')


plt.figure()
ax1= plt.axes()
ri = magR-magI
iz = magI - magZ
ax1.set_xlabel('r-i')
ax1.set_ylabel('i-z')
randomPlot(ri, iz)
plt.savefig('data/overlay/rivsiz.png')



plt.figure()
ax1= plt.axes()
iz = magI - magZ
zy = magZ - magY
ax1.set_xlabel('i-z')
ax1.set_ylabel('z-y')
randomPlot(iz, zy)
plt.savefig('data/overlay/izvszy.png')


plt.figure()
ax1= plt.axes()
zy = magZ - magY
yj = magY - jmags
ax1.set_ylabel('y-j')
randomPlot(zy, yj)
plt.savefig('data/overlay/zyvsyj.png')


plt.figure()
ax1= plt.axes()
yj = magY - jmags
jh = jmags - hmags
ax1.set_xlabel('y-j')
ax1.set_ylabel('j-h')
randomPlot(yj, jh)
plt.savefig('data/overlay/yjvsjh.png')

plt.figure()
ax1= plt.axes()
jh = jmags - hmags
hk = hmags - kmags
ax1.set_xlabel('j-h')
ax1.set_ylabel('h-k')
randomPlot(jh, hk)
plt.savefig('data/overlay/jhvshk.png')


plt.figure()
ax1= plt.axes()
hk = hmags - kmags
kchan1 = kmags - chan1mags
ax1.set_xlabel(r'$h-k$')
ax1.set_ylabel(r'$k-3.6\mu{}m$')
randomPlot(hk, kchan1)
plt.savefig('data/overlay/hkvsk36u.png')


plt.figure()
ax1= plt.axes()
kchan1 = kmags - chan1mags
chan1chan2 = chan1mags-chan2mags
ax1.set_xlabel(r'$k-3.6\mu{}m$')
ax1.set_ylabel(r'$3.6\mu{}m-4.5\mu{}m$')
randomPlot(kchan1, chan1chan2)
plt.savefig('data/overlay/kchan1vschan1chan2.png')


plt.figure()
ax1= plt.axes()
chan1chan2 = chan1mags-chan2mags
chan2chan3 = chan2mags-chan3mags
ax1.set_xlabel(r'$3.6-4.5\mu{}m$')
ax1.set_ylabel(r'$4.5\mu{}m-5.8\mu{}m$')
randomPlot(chan1chan2, chan2chan3)
plt.savefig('data/overlay/chan1chan2vschan2chan3.png')

plt.figure()
ax1= plt.axes()
chan2chan3 = chan2mags-chan3mags
chan3chan4 = chan3mags-chan4mags
ax1.set_xlabel(r'$4.5-5.8\mu{}m$')
ax1.set_ylabel(r'$5.8\mu{}m-8.0\mu{}m$')
randomPlot(chan2chan3, chan3chan4)
plt.savefig('data/overlay/chan2chan3vschan3chan4.png')

