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

badStarIndices = np.genfromtxt('data/deconv/starBadPredictionIndexProbability.txt')[:, 0].tolist()
badGalaxyIndices = np.genfromtxt('data/deconv/galaxyBadPredictionIndexProbability.txt')[:, 0].tolist()

magGBad = matchedCatBad['magG']
magRBad = matchedCatBad['magR']
magIBad = matchedCatBad['magI']
magZBad = matchedCatBad['magZ']
magYBad = matchedCatBad['magY']

ymagsBad = matchedCatBad['mag_y']
jmagsBad = matchedCatBad['mag_j']
hmagsBad = matchedCatBad['mag_h']
kmagsBad = matchedCatBad['mag_k']

chan1magsBad = matchedCatBad['mag_36']
chan2magsBad = matchedCatBad['mag_45']
chan3magsBad = matchedCatBad['mag_58']
chan4magsBad = matchedCatBad['mag_80']


def randomPlot(mag1, mag2, mag3, mag1Error, mag2Error, mag3Error, mag1bad, mag2bad, mag3bad, xaxistitle, yaxistitle, filename):
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
    

    starIndicesLocal = np.random.choice(starIndicesLocal, int(len(starIndicesLocal)*0.2), replace=False)
    galaxyIndicesLocal = np.random.choice(galaxyIndicesLocal, int(len(galaxyIndicesLocal) * 0.2), replace=False)

    print 'Number All Stars {0}'.format(len(starIndicesLocal))
    print 'Number All Galaxies {0}'.format(len(galaxyIndicesLocal))


    plt.scatter(magSource1[galaxyIndicesLocal], magSource2[galaxyIndicesLocal], marker='.', c='blue', \
            label = 'Galaxies')
    plt.scatter(magSource1[starIndicesLocal], magSource2[starIndicesLocal], marker='.', c='red', \
            label='Stars')
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles, labels)

    ax1.set_ylim([-2, 2])
    ax1.set_xlim([-2, 2])
    ax1.set_xlabel(xaxistitle, fontdict={'fontsize': 10})
    ax1.set_ylabel(yaxistitle, fontdict={'fontsize': 10})
    ax1.set_title('Distribution of Colors for All Sources')

    ax2 = plt.subplot(122)
    magSource1Bad = mag1bad - mag2bad
    magSource2Bad = mag2bad - mag3bad

    plt.scatter(magSource1Bad[badGalaxyIndices], magSource2Bad[badGalaxyIndices], marker='o', c='blue', \
            label = 'Galaxies')
    plt.scatter(magSource1Bad[badStarIndices], magSource2Bad[badStarIndices], marker='o', c='red', \
            label='Stars')

    ax2.set_ylim([-2, 2])
    ax2.set_xlim([-2, 2])
    ax2.set_xlabel(xaxistitle, fontdict={'fontsize': 10})
    ax2.set_ylabel(yaxistitle, fontdict={'fontsize': 10})
    ax2.set_title('Distribution of Colors for Misclassified Sources')
    handles, labels = ax2.get_legend_handles_labels()
    ax2.legend(handles, labels)


    # speed up performance
    fig.tight_layout()
    plt.savefig(filename)
    plt.close(fig)
    print 'done'

def randomPlotColors(color1, color2, mag1Error, mag2Error, mag3Error, mag4Error, color1Bad, color2Bad, xaxistitle, yaxistitle, filename):
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

    print 'Number All Stars {0}'.format(len(starIndicesLocal))
    print 'Number All Galaxies {0}'.format(len(galaxyIndicesLocal))

    plt.scatter(magSource1[galaxyIndicesLocal], magSource2[galaxyIndicesLocal], marker='o', c='blue', \
            label = 'Galaxies')
    plt.scatter(magSource1[starIndicesLocal], magSource2[starIndicesLocal], marker='o', c='red', \
            label='Stars')
    
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles, labels)

    ax1.set_ylim([-2, 2])
    ax1.set_xlim([-2, 2])
    ax1.set_xlabel(xaxistitle, fontdict={'fontsize': 10})
    ax1.set_ylabel(yaxistitle, fontdict={'fontsize': 10})
    ax1.set_title('Distribution of Colors for All Sources')

    ax2 = plt.subplot(122)
    magSource1Bad = color1Bad
    magSource2Bad = color2Bad

    print 'Number Bad Stars {0}'.format(len(badStarIndices))
    print 'Number Bad Galaxies {0}'.format(len(badGalaxyIndices))

    plt.scatter(magSource1Bad[badGalaxyIndices], magSource2Bad[badGalaxyIndices], c='blue', marker='o',\
            label = 'Galaxies')
    plt.scatter(magSource1Bad[badStarIndices], magSource2Bad[badStarIndices], c='red', marker='o',\
            label='Stars')

    ax2.set_ylim([-2, 2])
    ax2.set_xlim([-2, 2])
    ax2.set_xlabel(xaxistitle, fontdict={'fontsize': 10})
    ax2.set_ylabel(yaxistitle, fontdict={'fontsize': 10})
    ax2.set_title('Distribution of Colors for Misclassified Sources')
    handles, labels = ax2.get_legend_handles_labels()
    ax2.legend(handles, labels)


    # speed up performance
    fig.tight_layout()
    plt.savefig(filename)
    plt.close(fig)
    print 'done'


randomPlot(magGAll, magRAll, magIAll, magGErrorAll, magRErrorAll, magIErrorAll, magGBad, magRBad, magIBad,\
        'g-r', 'r-i', 'data/SideBySide/gri.png')
randomPlot(magRAll, magIAll, magZAll, magRErrorAll, magIErrorAll, magZErrorAll, magRBad, magIBad, magZBad, \
        'r-i', 'i-z', 'data/SideBySide/riz.png')
randomPlot(jmagsAll, hmagsAll, kmagsAll, matchedCatAll['mag_j_error'], matchedCatAll['mag_h_error'], matchedCatAll['mag_k_error'],\
        jmagsBad, hmagsBad, kmagsBad, \
        'j-h', 'h-k', 'data/SideBySide/jhk.png')

randomPlotColors(hmagsAll-kmagsAll, chan1magsAll-chan2magsAll, matchedCatAll['mag_h_error'], matchedCatAll['mag_k_error'], \
                matchedCatAll['mag_36error'], matchedCatAll['mag_45error'], hmagsBad - kmagsBad, chan1magsBad-chan2magsBad, \
                r'$3.6\mu{}m - 4.5\mu{}m$', 'h-k', 'data/SideBySide/hkchan1chan2.png')
