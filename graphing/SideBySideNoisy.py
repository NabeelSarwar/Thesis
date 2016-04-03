import matplotlib
import gc
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.markers as markers
plt.ioff()

import numpy as np
import pyfits
matchedCatAll = pyfits.getdata('MatchHell2.fits')
#shortcut and
sand = np.logical_and
markerS = markers.MarkerStyle(marker='.')
goodIndices = matchedCatAll['cmodel_flux_g'] > 0
goodIndices = sand(goodIndices, matchedCatAll['cmodel_flux_r'] > 0)
goodIndices = sand(goodIndices, matchedCatAll['cmodel_flux_i'] > 0)
goodIndices = sand(goodIndices, matchedCatAll['cmodel_flux_z'] > 0)
goodIndices = sand(goodIndices, matchedCatAll['cmodel_flux_y'] > 0)
goodIndices = sand(goodIndices, matchedCatAll['cmodel_flux_err_g'] != 0)
goodIndices = sand(goodIndices, matchedCatAll['cmodel_flux_err_r'] != 0)
goodIndices = sand(goodIndices, matchedCatAll['cmodel_flux_err_i'] != 0)
goodIndices = sand(goodIndices, matchedCatAll['cmodel_flux_err_z'] != 0)
goodIndices = sand(goodIndices, matchedCatAll['cmodel_flux_err_y'] != 0)

matchedCatAll = pyfits.BinTableHDU(matchedCatAll[goodIndices]).data

starindicesAll = matchedCatAll['mu_class'] == 2
galaxyindicesAll = matchedCatAll['mu_class'] == 1

magGAll = -2.5*np.log10(matchedCatAll['cmodel_flux_g']/matchedCatAll['flux_zeromag_g'])
magRAll = -2.5*np.log10(matchedCatAll['cmodel_flux_r']/matchedCatAll['flux_zeromag_r'])
magIAll = -2.5*np.log10(matchedCatAll['cmodel_flux_i']/matchedCatAll['flux_zeromag_i'])
magZAll = -2.5*np.log10(matchedCatAll['cmodel_flux_z']/matchedCatAll['flux_zeromag_z'])
magYAll = -2.5*np.log10(matchedCatAll['cmodel_flux_y']/matchedCatAll['flux_zeromag_y'])

magGErrorAll = np.abs(1.08574 / matchedCatAll['cmodel_flux_g'] * matchedCatAll['cmodel_flux_err_g'])
magRErrorAll = np.abs(1.08574 / matchedCatAll['cmodel_flux_r'] * matchedCatAll['cmodel_flux_err_r'])
magIErrorAll = np.abs(1.08574 / matchedCatAll['cmodel_flux_i'] * matchedCatAll['cmodel_flux_err_i'])
magZErrorAll = np.abs(1.08574 / matchedCatAll['cmodel_flux_z'] * matchedCatAll['cmodel_flux_err_z'])
magYErrorAll = np.abs(1.08574 / matchedCatAll['cmodel_flux_y'] * matchedCatAll['cmodel_flux_err_y'])

ymagsAll = matchedCatAll['mag_y']
jmagsAll = matchedCatAll['mag_j']
hmagsAll = matchedCatAll['mag_h']
kmagsAll = matchedCatAll['mag_k']

chan1magsAll = matchedCatAll['mag_36']
chan2magsAll = matchedCatAll['mag_45']
chan3magsAll = matchedCatAll['mag_58']
chan4magsAll = matchedCatAll['mag_80']

# now let's get the bad predictions and such
matchedCatBad = pyfits.getdata('MatchHellDeconv2.fits')

# yeah this is bad, but i am in a rush
badStarIndicesGlob = np.genfromtxt('data/deconvnoisy/starBadPredictionIndexProbability.txt')[:, 0].tolist()
badGalaxyIndicesGlob = np.genfromtxt('data/deconvnoisy/galaxyBadPredictionIndexProbability.txt')[:, 0].tolist()

magGBad = matchedCatBad['magG']
magRBad = matchedCatBad['magR']
magIBad = matchedCatBad['magI']
magZBad = matchedCatBad['magZ']
magYBad = matchedCatBad['magY']

cleanDeconvG = matchedCatBad['magGError'] < 0.2
cleanDeconvR = matchedCatBad['magRError'] < 0.2
cleanDeconvI = matchedCatBad['magIError'] < 0.2
cleanDeconvZ = matchedCatBad['magZError'] < 0.2

ymagsBad = matchedCatBad['mag_y']
jmagsBad = matchedCatBad['mag_j']
hmagsBad = matchedCatBad['mag_h']
kmagsBad = matchedCatBad['mag_k']

cleanDeconvJ = matchedCatBad['mag_j_error'] < 0.2
cleanDeconvH = matchedCatBad['mag_h_error'] < 0.2
cleanDeconvK = matchedCatBad['mag_k_error'] < 0.2

chan1magsBad = matchedCatBad['mag_36']
chan2magsBad = matchedCatBad['mag_45']
chan3magsBad = matchedCatBad['mag_58']
chan4magsBad = matchedCatBad['mag_80']

cleanDeconvChan1 = matchedCatBad['mag_36error'] < 0.2
cleanDeconvChan2 = matchedCatBad['mag_45error'] < 0.2
cleanDeconvChan3 = matchedCatBad['mag_58error'] < 0.2
cleanDeconvChan4 = matchedCatBad['mag_80error'] < 0.2


def randomPlot(mag1, mag2, mag3, mag1Error, mag2Error, mag3Error, mag1bad, mag2bad, mag3bad, \
        cleanIndicesBad1, cleanIndicesBad2, cleanIndicesBad3, xaxistitle, yaxistitle, filename):
    fig = plt.figure(figsize=(10, 5))

    ax1 = plt.subplot(121)

    magSource1 = mag1 - mag2
    magSource2 = mag2 - mag3
    good12 = sand(mag1Error < 0.2, mag2Error < 0.2)
    good23 = sand(mag2Error < 0.2, mag3Error <0.2)
    goodIndices = sand(good12, good23)

    starIndicesLocal = sand(starindicesAll, goodIndices)
    galaxyIndicesLocal = sand(galaxyindicesAll, goodIndices)

    # because I sample indices
    starIndicesLocal = np.where(starIndicesLocal==True)[0]
    galaxyIndicesLocal = np.where(galaxyIndicesLocal==True)[0]
    

    # take a random sample of our data since there is so much
    starIndicesLocal = np.random.choice(starIndicesLocal, int(len(starIndicesLocal)*0.2), replace=False)
    galaxyIndicesLocal = np.random.choice(galaxyIndicesLocal, int(len(galaxyIndicesLocal) * 0.2), replace=False)

    plt.scatter(magSource1[galaxyIndicesLocal], magSource2[galaxyIndicesLocal], marker='.', c='blue', \
            label = 'Galaxies')
    plt.scatter(magSource1[starIndicesLocal], magSource2[starIndicesLocal], marker='.', c='red', \
            label='Stars')
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles, labels, loc='bottom right')

    ax1.set_ylim([-2, 2])
    ax1.set_xlim([-2, 2])
    ax1.set_xlabel(xaxistitle, fontdict={'fontsize': 10})
    ax1.set_ylabel(yaxistitle, fontdict={'fontsize': 10})
    ax1.set_title('Distribution of Colors for All Sources')

    ax2 = plt.subplot(122)
    magSource1Bad = mag1bad - mag2bad
    magSource2Bad = mag2bad - mag3bad

    # bad galaxy indices and bad star indices are index counts, convert them to mask
    # mag source 1 bad should have the right number of entries

    badStarIndices = [True if i in badStarIndicesGlob else False for i in range(len(magSource1Bad))]
    badGalaxyIndices = [True if i in badGalaxyIndicesGlob else False for i in range(len(magSource1Bad))]
    cleanSource1Bad = np.logical_and(cleanIndicesBad1, cleanIndicesBad2)
    cleanSource2Bad = np.logical_and(cleanIndicesBad2, cleanIndicesBad3)

    plt.scatter(magSource1Bad[sand(badGalaxyIndices, cleanSource1Bad)], \
                magSource2Bad[sand(badGalaxyIndices, cleanSource2Bad)], \
                marker='.', c='blue', label = 'Galaxies')
    plt.scatter(magSource1Bad[sand(badStarIndices, cleanSource1Bad)], \
                magSource2Bad[sand(badStarIndices, cleanSource2Bad)], \
                marker='.', c='red', label='Stars')

    # plot the noisy data
    plt.scatter(magSource1Bad[sand(badGalaxyIndices, np.logical_not(cleanSource1Bad))], \
                magSource2Bad[sand(badGalaxyIndices, np.logical_not(cleanSource2Bad))], \
                marker='s', c='blue', label = 'Noisy Galaxies')
    plt.scatter(magSource1Bad[sand(badStarIndices, np.logical_not(cleanSource1Bad))], \
                magSource2Bad[sand(badStarIndices, np.logical_not(cleanSource2Bad))], \
                marker='s', c='red', label='Noisy Stars')


    ax2.set_ylim([-2, 2])
    ax2.set_xlim([-2, 2])
    ax2.set_xlabel(xaxistitle, fontdict={'fontsize': 10})
    ax2.set_ylabel(yaxistitle, fontdict={'fontsize': 10})
    ax2.set_title('Distribution of Colors for Misclassified Sources')
    handles, labels = ax2.get_legend_handles_labels()
    ax2.legend(handles, labels, loc='bottom right')


    # speed up performance
    fig.tight_layout()
    plt.savefig(filename)
    plt.close(fig)
    print 'done'

def randomPlotColors(color1, color2, mag1Error, mag2Error, mag3Error, mag4Error, color1Bad, color2Bad, \
        cleanIndicesBad1, cleanIndicesBad2, cleanIndicesBad3, cleanIndicesBad4, \
        xaxistitle, yaxistitle, filename):
    fig = plt.figure(figsize=(10,5))

    ax1 = plt.subplot(121)

    magSource1 = color1
    magSource2 = color2
    good12 = sand(mag1Error < 0.2, mag2Error < 0.2)
    good23 = sand(mag3Error < 0.2, mag4Error <0.2)
    goodIndices = sand(good12, good23)

    starIndicesLocal = sand(starindicesAll, goodIndices)
    galaxyIndicesLocal = sand(galaxyindicesAll, goodIndices)


    # because I sample indices
    starIndicesLocal = np.where(starIndicesLocal==True)[0]
    galaxyIndicesLocal = np.where(galaxyIndicesLocal==True)[0]
 
    starIndicesLocal = np.random.choice(starIndicesLocal, int(len(starIndicesLocal)*0.2), replace=False)
    galaxyIndicesLocal = np.random.choice(galaxyIndicesLocal, int(len(galaxyIndicesLocal) * 0.2), replace=False)


    plt.scatter(magSource1[galaxyIndicesLocal], magSource2[galaxyIndicesLocal], marker='.', c='blue', \
            label = 'Galaxies')
    plt.scatter(magSource1[starIndicesLocal], magSource2[starIndicesLocal], marker='.', c='red', \
            label='Stars')
    
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles, labels, loc='bottom right')

    ax1.set_ylim([-2, 2])
    ax1.set_xlim([-2, 2])
    ax1.set_xlabel(xaxistitle, fontdict={'fontsize': 10})
    ax1.set_ylabel(yaxistitle, fontdict={'fontsize': 10})
    ax1.set_title('Distribution of Colors for All Sources')

    ax2 = plt.subplot(122)
    magSource1Bad = color1Bad
    magSource2Bad = color2Bad

    badStarIndices = [True if i in badStarIndicesGlob else False for i in range(len(magSource1Bad))]
    badGalaxyIndices = [True if i in badGalaxyIndicesGlob else False for i in range(len(magSource1Bad))]


    # separate the noisy stars and galaxies from the non noisy ones (all from improper classifications)
    indices1 = sand(badGalaxyIndices, sand(cleanIndicesBad1, cleanIndicesBad2))
    indices2 = sand(badGalaxyIndices, sand(cleanIndicesBad3, cleanIndicesBad4))
    indices3 = sand(badStarIndices, sand(cleanIndicesBad1, cleanIndicesBad2))
    indices4 = sand(badStarIndices, sand(cleanIndicesBad3, cleanIndicesBad4))

    plt.scatter(magSource1Bad[indices1], magSource2Bad[indices2], c='blue', marker='.',\
            label = 'Galaxies')
    plt.scatter(magSource1Bad[indices3], magSource2Bad[indices4], c='red', marker='.',\
            label='Stars')

    indices1 = sand(badGalaxyIndices, np.logical_not(sand(cleanIndicesBad1, cleanIndicesBad2)))
    indices2 = sand(badGalaxyIndices, np.logical_not(sand(cleanIndicesBad3, cleanIndicesBad4)))
    indices3 = sand(badStarIndices, np.logical_not(sand(cleanIndicesBad1, cleanIndicesBad2)))
    indices4 = sand(badStarIndices, np.logical_not(sand(cleanIndicesBad3, cleanIndicesBad4)))

    plt.scatter(magSource1Bad[indices1], magSource2Bad[indices2], c='blue', marker='s',\
            label = 'Noisy Galaxies')
    plt.scatter(magSource1Bad[indices3], magSource2Bad[indices4], c='red', marker='s',\
            label='Noisy Stars')

    ax2.set_ylim([-2, 2])
    ax2.set_xlim([-2, 2])
    ax2.set_xlabel(xaxistitle, fontdict={'fontsize': 10})
    ax2.set_ylabel(yaxistitle, fontdict={'fontsize': 10})
    ax2.set_title('Distribution of Colors for Misclassified Sources')
    handles, labels = ax2.get_legend_handles_labels()
    ax2.legend(handles, labels, loc='bottom right')


    # speed up performance
    fig.tight_layout()
    plt.savefig(filename)
    plt.close(fig)
    print 'done'


randomPlot(magGAll, magRAll, magIAll, magGErrorAll, magRErrorAll, magIErrorAll, magGBad, magRBad, magIBad,\
        cleanDeconvG, cleanDeconvR, cleanDeconvI, \
        'g-r', 'r-i', 'data/SideBySide/gri.png')
randomPlot(magRAll, magIAll, magZAll, magRErrorAll, magIErrorAll, magZErrorAll, magRBad, magIBad, magZBad, \
        cleanDeconvR, cleanDeconvI, cleanDeconvZ, \
        'r-i', 'i-z', 'data/SideBySide/riz.png')
randomPlot(jmagsAll, hmagsAll, kmagsAll, matchedCatAll['mag_j_error'], matchedCatAll['mag_h_error'], matchedCatAll['mag_k_error'],\
        jmagsBad, hmagsBad, kmagsBad, \
        cleanDeconvJ, cleanDeconvH, cleanDeconvK, 
        'j-h', 'h-k', 'data/SideBySide/jhk.png')

randomPlotColors(hmagsAll-kmagsAll, chan1magsAll-chan2magsAll, matchedCatAll['mag_h_error'], matchedCatAll['mag_k_error'], \
                matchedCatAll['mag_36error'], matchedCatAll['mag_45error'], hmagsBad - kmagsBad, chan1magsBad-chan2magsBad, \
                cleanDeconvH, cleanDeconvK, cleanDeconvChan1, cleanDeconvChan2, \
                r'$3.6\mu{}m - 4.5\mu{}m$', 'h-k', 'data/SideBySide/hkchan1chan2.png')
