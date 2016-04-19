from copy import deepcopy
import pyfits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, Column
from astropy.io import ascii

# get magnitudes from match deep

matchdeepcatalog = pyfits.getdata('matchDeepCoaddMeas-137520151126CosmosGRIZY.fits')

magG = -2.5*np.log10(matchdeepcatalog['cmodel_flux_g']/matchdeepcatalog['flux_zeromag_g'])
magR = -2.5*np.log10(matchdeepcatalog['cmodel_flux_r']/matchdeepcatalog['flux_zeromag_r'])
magI = -2.5*np.log10(matchdeepcatalog['cmodel_flux_i']/matchdeepcatalog['flux_zeromag_i'])
magZ = -2.5*np.log10(matchdeepcatalog['cmodel_flux_z']/matchdeepcatalog['flux_zeromag_z'])
magY = -2.5*np.log10(matchdeepcatalog['cmodel_flux_y']/matchdeepcatalog['flux_zeromag_y'])

starIndices = np.where(matchdeepcatalog['mu_class']==2)[0]
galaxyIndices = np.where(matchdeepcatalog['mu_class']==1)[0]

def PlotHistogramMagnitudes(magnitudes, bins, fileName, title):
    binCounts = np.zeros(len(bins))
    bins2 = np.append(bins, [np.infty])

    magnitudes = magnitudes[np.logical_not(np.isnan(magnitudes))]
    for i in range(len(bins2)-1):
        binCounts[i] = np.sum(np.logical_and(magnitudes >= bins2[i], magnitudes < bins2[i+1]))

    fig = plt.figure()

    plt.title(title)
    plt.xlabel('Magnitude (AB)')
    plt.ylabel(r'$\log_{10}(\Delta{}N)$')
    binCounts = np.log10(binCounts)
    width = 1.0
    plt.xticks(bins+width, bins)
    plt.bar(bins, binCounts, width)
    plt.savefig(fileName)

bins = np.arange(30)

#HSC
#Stars only
PlotHistogramMagnitudes(magR[starIndices], bins, 'data/logcountsall/countsStarsRBand.png', r'Logarithm of Approximate $\frac{dN_{star}}{dm}$ in R Band (AB)')
PlotHistogramMagnitudes(magI[starIndices], bins, 'data/logcountsall/countsStarsIBand.png', r'Logarithm of Approximate $\frac{dN_{star}}{dm}$ in I Band (AB)')
PlotHistogramMagnitudes(magZ[starIndices], bins, 'data/logcountsall/countsStarsZBand.png', r'Logarithm of Approximate $\frac{dN_{star}}{dm}$ in Z Band (AB)')
PlotHistogramMagnitudes(magY[starIndices], bins, 'data/logcountsall/countsStarsYBand.png', r'Logarithm of Approximate $\frac{dN_{star}}{dm}$ in Y Band (AB)')



#Galaxies
PlotHistogramMagnitudes(magR[galaxyIndices], bins, 'data/logcountsall/countsGalaxiesRBand.png', r'Logarithm of Approximate $\frac{dN_{galaxy}}{dm}$ in R Band (AB)')
PlotHistogramMagnitudes(magI[galaxyIndices], bins, 'data/logcountsall/countsGalaxiesIBand.png', r'Logarithm of Approximate $\frac{dN_{galaxy}}{dm}$ in I Band (AB)')
PlotHistogramMagnitudes(magZ[galaxyIndices], bins, 'data/logcountsall/countsGalaxiesZBand.png', r'Logarithm of Approximate $\frac{dN_{galaxy}}{dm}$ in Z Band (AB)')
PlotHistogramMagnitudes(magY[galaxyIndices], bins, 'data/logcountsall/countsStarsYBand.png', r'Logarithm of Approximate $\frac{dN_{star}}{dm}$ in Y Band (AB)')



#Ultravista
data = pyfits.getdata('UltraVistaHSCHST.fits')
starIndices = np.where(data['mu_class']==2)[0]
galaxyIndices = np.where(data['mu_class']==1)[0]

magj = data['mag_j']
magh = data['mag_h']
magk = data['mag_k']

# Stars
PlotHistogramMagnitudes(magj[starIndices], bins, 'data/logcountsall/countsStarsJBand.png', r'Logarithm of Approximate $\frac{dN}{dm}$ in J Band (AB)')
PlotHistogramMagnitudes(magh[starIndices], bins, 'data/logcountsall/countsStarsHBand.png', r'Logarithm of Approximate $\frac{dN}{dm}$ in H Band (AB)')
PlotHistogramMagnitudes(magk[starIndices], bins, 'data/logcountsall/countsStarsKBand.png', r'Logarithm of Approximate $\frac{dN}{dm}$ in K Band (AB)')

# Galaxies
PlotHistogramMagnitudes(magj[galaxyIndices], bins, 'data/logcountsall/countsGalaxiesJBand.png', r'Logarithm of Approximate $\frac{dN}{dm}$ in J Band (AB)')
PlotHistogramMagnitudes(magh[galaxyIndices], bins, 'data/logcountsall/countsGalaxiesHBand.png', r'Logarithm of Approximate $\frac{dN}{dm}$ in H Band (AB)')
PlotHistogramMagnitudes(magk[galaxyIndices], bins, 'data/logcountsall/countsGalaxiesKBand.png', r'Logarithm of Approximate $\frac{dN}{dm}$ in K Band (AB)')


#Spitzer
data = pyfits.getdata('SpitzerHSCHST.fits')
starIndices = np.where(data['mu_class']==2)[0]
galaxyIndices = np.where(data['mu_class']==1)[0]

mag36 = data['mag_36'] 
mag45 = data['mag_45']

# Stars
PlotHistogramMagnitudes(mag36[starIndices], bins, 'data/logcountsall/countsStars36MicronBand.png', r'Logarithm of Approximate $\frac{dN}{dm}$ in 3.6 Band (AB)')
PlotHistogramMagnitudes(mag45[starIndices], bins, 'data/logcountsall/countsStars45MicronBand.png', r'Logarithm of Approximate $\frac{dN}{dm}$ in 4.5 Band (AB)')

# Galaxies
PlotHistogramMagnitudes(mag36[galaxyIndices], bins, 'data/logcountsall/countsGalaxies36MicronBand.png', r'Logarithm of Approximate $\frac{dN}{dm}$ in 3.6 Band (AB)')
PlotHistogramMagnitudes(mag45[galaxyIndices], bins, 'data/logcountsall/countsGalaxies45MicronBand.png', r'Logarithm of Approximate $\frac{dN}{dm}$ in 4.5 Band (AB)')
