import matplotlib
import gc
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.markers as markers
plt.ioff()

import numpy as np
import pyfits
matchedCat = pyfits.getdata('MatchHell.fits')
#shortcut and
sand = np.logical_and
markerS = markers.MarkerStyle(marker='.')
goodIndices = matchedCat['cmodel_flux_g'] > 0
goodIndices = sand(goodIndices, matchedCat['cmodel_flux_r'] > 0)
goodIndices = sand(goodIndices, matchedCat['cmodel_flux_i'] > 0)
goodIndices = sand(goodIndices, matchedCat['cmodel_flux_z'] > 0)
goodIndices = sand(goodIndices, matchedCat['cmodel_flux_y'] > 0)
goodIndices = sand(goodIndices, matchedCat['cmodel_flux_err_g'] != 0)
goodIndices = sand(goodIndices, matchedCat['cmodel_flux_err_r'] != 0)
goodIndices = sand(goodIndices, matchedCat['cmodel_flux_err_i'] != 0)
goodIndices = sand(goodIndices, matchedCat['cmodel_flux_err_z'] != 0)
goodIndices = sand(goodIndices, matchedCat['cmodel_flux_err_y'] != 0)
matchedCat = matchedCat[goodIndices]
starindices = np.where((matchedCat['mu_class'] == 2)==True)[0]
galaxyindices = np.where((matchedCat['mu_class'] == 1)==True)[0]
print starindices
print galaxyindices
starGalHash = {}
for e in starindices:
    starGalHash[e] = 'star'

for e in galaxyindices:
    starGalHash[e] = 'galaxy'

magG = -2.5*np.log10(matchedCat['cmodel_flux_g']/matchedCat['flux_zeromag_g'])
magR = -2.5*np.log10(matchedCat['cmodel_flux_r']/matchedCat['flux_zeromag_r'])
magI = -2.5*np.log10(matchedCat['cmodel_flux_i']/matchedCat['flux_zeromag_i'])
magZ = -2.5*np.log10(matchedCat['cmodel_flux_z']/matchedCat['flux_zeromag_z'])
magY = -2.5*np.log10(matchedCat['cmodel_flux_y']/matchedCat['flux_zeromag_y'])

magGError = np.abs(1.08574 / matchedCat['cmodel_flux_g'] * matchedCat['cmodel_flux_err_g'])
magRError = np.abs(1.08574 / matchedCat['cmodel_flux_r'] * matchedCat['cmodel_flux_err_r'])
magIError = np.abs(1.08574 / matchedCat['cmodel_flux_i'] * matchedCat['cmodel_flux_err_i'])
magZError = np.abs(1.08574 / matchedCat['cmodel_flux_z'] * matchedCat['cmodel_flux_err_z'])
magYError = np.abs(1.08574 / matchedCat['cmodel_flux_y'] * matchedCat['cmodel_flux_err_y'])

ymags = matchedCat['mag_y']
jmags = matchedCat['mag_j']
hmags = matchedCat['mag_h']
kmags = matchedCat['mag_k']

chan1mags = matchedCat['mag_36']
chan2mags = matchedCat['mag_45']
chan3mags = matchedCat['mag_58']
chan4mags = matchedCat['mag_80']

def randomPlot(magSource1, magSource2, goodIndices):
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
goodGR = sand(magGError < 0.2, magRError < 0.2)
goodRI = sand(magRError < 0.2, magIError < 0.2)
good = sand(goodGR, goodRI)
ax1.set_xlabel('g-r')
ax1.set_ylabel('r-i')
randomPlot(gr, ri, good)
#plt.scatter(gr[sand(starindices, good)], ri[sand(starindices, good)], marker='+', c='red', label='star')
#plt.scatter(gr[sand(galaxyindices, good)], ri[sand(galaxyindices, good)], marker='+', c='blue', label='galaxy')
plt.savefig('colorplots/grvsri.png')


plt.figure()
ax1= plt.axes()
ri = magR-magI
iz = magI - magZ
goodIZ = sand(magIError < 0.2, magZError < 0.2)
good = sand(goodRI, goodIZ)
ax1.set_xlabel('r-i')
ax1.set_ylabel('i-z')
randomPlot(ri, iz, good)
#plt.scatter(ri[sand(starindices, good)], iz[sand(starindices, good)], marker='+', c='red', label='star')
#plt.scatter(ri[sand(galaxyindices, good)], iz[sand(galaxyindices, good)], marker='+', c='blue', label='galaxy')
plt.savefig('colorplots/rivsiz.png')



plt.figure()
ax1= plt.axes()
iz = magI - magZ
zy = magZ - magY
goodZY = sand(magZError < 0.2, magYError < 0.2)
good = sand(goodIZ, goodZY)
ax1.set_xlabel('i-z')
ax1.set_ylabel('z-y')
randomPlot(iz, zy, good)
#plt.scatter(iz[sand(starindices, good)], zy[sand(starindices, good)], marker='+', c='red', label='star')
#plt.scatter(iz[sand(galaxyindices, good)], zy[sand(galaxyindices, good)], marker='+', c='blue', label='galaxy')
plt.savefig('colorplots/izvszy.png')


plt.figure()
ax1= plt.axes()
zy = magZ - magY
yj = magY - jmags
goodYJ = sand(magYError < 0.2, matchedCat['mag_j_error'] < 0.2)
good = sand(goodZY, goodYJ)
ax1.set_xlabel('z-y')
ax1.set_ylabel('y-j')
randomPlot(zy, yj, good)
#plt.scatter(zy[sand(starindices, good)], yj[sand(starindices, good)], marker='+', c='red', label='star')
#plt.scatter(zy[sand(galaxyindices, good)], yj[sand(galaxyindices, good)], marker='+', c='blue', label='galaxy')
plt.savefig('colorplots/zyvsyj.png')


plt.figure()
ax1= plt.axes()
yj = magY - jmags
jh = jmags - hmags
goodJH = sand(matchedCat['mag_j_error'] < 0.2, matchedCat['mag_h_error'] < 0.2)
good = sand(goodYJ, goodJH)
ax1.set_xlabel('y-j')
ax1.set_ylabel('j-h')
randomPlot(yj, jh, good)
#plt.scatter(yj[sand(starindices, good)], jh[sand(starindices, good)], marker='+', c='red', label='star')
#plt.scatter(yj[sand(galaxyindices, good)], jh[sand(galaxyindices, good)], marker='+', c='blue', label='galaxy')
plt.savefig('colorplots/yjvsjh.png')

plt.figure()
ax1= plt.axes()
jh = jmags - hmags
hk = hmags - kmags
goodHK = sand(matchedCat['mag_h_error'] < 0.2, matchedCat['mag_k_error'] < 0.2)
good = sand(goodJH, goodHK)
ax1.set_xlabel('j-h')
ax1.set_ylabel('h-k')
randomPlot(jh, hk, good)
#plt.scatter(jh[sand(starindices, good)], hk[sand(starindices, good)], marker='+', c='red', label='star')
#plt.scatter(jh[sand(galaxyindices, good)], hk[sand(galaxyindices, good)], marker='+', c='blue', label='galaxy')
plt.savefig('colorplots/jhvshk.png')


plt.figure()
ax1= plt.axes()
hk = hmags - kmags
kchan1 = kmags - chan1mags
goodKChan1 = sand(matchedCat['mag_k_error'] < 0.2, matchedCat['mag_36error'] < 0.2)
good = sand(goodHK, goodKChan1)
ax1.set_xlabel(r'$h-k$')
ax1.set_ylabel(r'$k-3.6\mu{}m$')
randomPlot(hk, kchan1, good)
#plt.scatter(hk[sand(starindices, good)], kchan1[sand(starindices, good)], marker='+', c='red', label='star')
#plt.scatter(hk[sand(galaxyindices, good)], kchan1[sand(galaxyindices, good)], marker='+', c='blue', label='galaxy')
plt.savefig('colorplots/hkvsk36u.png')


plt.figure()
ax1= plt.axes()
kchan1 = kmags - chan1mags
chan1chan2 = chan1mags-chan2mags
goodChan1Chan2 = sand(matchedCat['mag_36error'] < 0.2, matchedCat['mag_45error'] < 0.2)
good = sand(goodKChan1, goodChan1Chan2)
ax1.set_xlabel(r'$k-3.6\mu{}m$')
ax1.set_ylabel(r'$3.6\mu{}m-4.5\mu{}m$')
randomPlot(kchan1, chan1chan2, good)
#plt.scatter(kchan1[sand(starindices, good)], chan1chan2[sand(starindices, good)], marker='+', c='red', label='star')
#plt.scatter(kchan1[sand(galaxyindices, good)], chan1chan2[sand(galaxyindices, good)], marker='+', c='blue', label='galaxy')
plt.savefig('colorplots/kchan1vschan1chan2.png')


plt.figure()
ax1= plt.axes()
chan1chan2 = chan1mags-chan2mags
chan2chan3 = chan2mags-chan3mags
goodChan2Chan3 = sand(matchedCat['mag_45error'] <0.2, matchedCat['mag_58error'] < 0.2)
good = sand(goodChan1Chan2, goodChan2Chan3)
ax1.set_xlabel(r'$3.6-4.5\mu{}m$')
ax1.set_ylabel(r'$4.5\mu{}m-5.8\mu{}m$')
randomPlot(chan1chan2, chan2chan3, good)
#plt.scatter(chan1chan2[sand(starindices, good)], chan2chan3[sand(starindices, good)], marker='+', c='red', label='star')
#plt.scatter(chan1chan2[sand(galaxyindices, good)], chan2chan3[sand(galaxyindices, good)], marker='+', c='blue', label='galaxy')
plt.savefig('colorplots/chan1chan2vschan2chan3.png')

plt.figure()
ax1= plt.axes()
chan2chan3 = chan2mags-chan3mags
chan3chan4 = chan3mags-chan4mags
goodChan3Chan4 = sand(matchedCat['mag_58error'] < 0.2, matchedCat['mag_80error'] < 0.2)
good = sand(goodChan2Chan3, goodChan3Chan4)
ax1.set_xlabel(r'$4.5-5.8\mu{}m$')
ax1.set_ylabel(r'$5.8\mu{}m-8.0\mu{}m$')
randomPlot(chan2chan3, chan3chan4, good)
#plt.scatter(chan2chan3[sand(starindices, good)], chan3chan4[sand(starindices, good)], marker='+', c='red', label='star')
#plt.scatter(chan2chan3[sand(galaxyindices, good)], chan3chan4[sand(galaxyindices, good)], marker='+', c='blue', label='galaxy')
plt.savefig('colorplots/chan2chan3vschan3chan4.png')
