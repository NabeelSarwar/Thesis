import astrML.density_estimation.XDGMM as XDGMM
import numpy as np
import scipy as sci
import matplotlib.pyplot as plt
import pyfits


ncomponents = 10
n_iter = 100
tol = 10**(-5)
clf = XDGMM.XDGMM(ncomponents, n_iter, tol)

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

#use this matchedCat for the proper compilation of matches
matchedCat = matchedCat[goodIndices]
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

# establish the indices with proper colors, and then get this
goodIndices = magRError < 0
goodIndices = sand(goodIndices, magIError < 0.2)
goodIndices = sand(goodIndices, magZError < 0.2)
goodIndices = sand(goodIndices, matchedCat['mag_j_error'] < 0.2)
goodIndices = sand(goodIndices, matchedCat['mag_h_error'] < 0.2)
goodIndices = sand(goodIndices, matchedCat['mag_k_error'] < 0.2)
goodIndices = sand(goodIndices, matchedCat['mag_36error'] < 0.2)
goodIndices = sand(goodIndices, matchedCat['mag_45error'] < 0.2)


# now  we have the proper errors
matchedCat = matchedCat[goodIndices]

ymags = matchedCat['mag_y']
jmags = matchedCat['mag_j']
hmags = matchedCat['mag_h']
kmags = matchedCat['mag_k']

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

chan1mags = matchedCat['mag_36']
chan2mags = matchedCat['mag_45']
chan3mags = matchedCat['mag_58']
chan4mags = matchedCat['mag_80']

starindices = np.where((matchedCat['mu_class'] == 2)==True)[0]
galaxyindices = np.where((matchedCat['mu_class'] == 1)==True)[0]
print starindices
print galaxyindices
starGalHash = {}
for e in starindices:
    starGalHash[e] = 'star'

for e in galaxyindices:
    starGalHash[e] = 'galaxy'


results = matchedCat['mu_class']
ids = matchedCat['id']

magnitudeMatrix = np.matrix([magR, magI, magZ, jmags, hmags, kmags, chan1mags, chan2mags], dtype=np.dtype('float64')).transpose()

# form the matrix of the good COLORS and the good Color errs
ricolor = magR-magI
izcolor = magI - magZ
jhcolor = jmags-hmags
hkcolor = hmags-kmags
kchan1color = kmags-chan1mags
chan1chan2color = chan1mags-chan2mags

def buildNoiseMatrix(catalog, iD):
    pass
