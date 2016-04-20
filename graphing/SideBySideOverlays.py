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

# this is needed for proper indexing because reading this way will give floats
deconvStarIndices = np.genfromtxt('data/deconv/starBadPredictionIndexProbability.txt')[:, 0].tolist()
deconvGalaxyIndices = np.genfromtxt('data/deconv/galaxyBadPredictionIndexProbability.txt')[:, 0].tolist()

magGDeconv = matchedCatBad['magG']
magRDeconv = matchedCatBad['magR']
magIDeconv = matchedCatBad['magI']
magZDeconv = matchedCatBad['magZ']
magYDeconv = matchedCatBad['magY']


ymagsDeconv = matchedCatBad['mag_y']
jmagsDeconv = matchedCatBad['mag_j']
hmagsDeconv = matchedCatBad['mag_h']
kmagsDeconv = matchedCatBad['mag_k']


chan1magsDeconv = matchedCatBad['mag_36']
chan2magsDeconv = matchedCatBad['mag_45']
chan3magsDeconv = matchedCatBad['mag_58']
chan4magsDeconv = matchedCatBad['mag_80']


# get the data from SVM and regression
matchedCat = pyfits.getdata('MatchHellErrorCut2.fits')
#shortcut and
sand = np.logical_and

def ErrorCutIDsToMask(idArray):
    indices = []
    for i in range(len(matchedCat['id'])):
        if matchedCat['id'][i] in idArray:
            indices.append(i)
    mask = np.array([True if i in indices else False for i in range(len(matchedCat['id']))])
    return mask

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
matchedCat = pyfits.BinTableHDU(matchedCat).data

idsStarLog = np.genfromtxt('data/logistic/badStarIDs.txt')[:, 0]
idsGalaxyLog = np.genfromtxt('data/logistic/badGalaxyIDs.txt')[:, 0]
idsStarLog = np.array([int(i) for i in idsStarLog])
idsGalaxyLog = np.array([int(i) for i in idsGalaxyLog])

idsStarSVML = np.genfromtxt('data/svm/badStarIDslinear.txt')[:, 0]
idsGalaxySVML = np.genfromtxt('data/svm/badGalaxyIDslinear.txt')[:, 0]
idsStarSVML = np.array([int(i) for i in idsStarSVML])
idsGalaxySVML = np.array([int(i) for i in idsGalaxySVML])

idsStarSVMR = np.genfromtxt('data/svm/badStarIDsRBF.txt')[:, 0]
idsGalaxySVMR = np.genfromtxt('data/svm/badGalaxyIDsRBF.txt')[:, 0]
idsStarSVMR = np.array([int(i) for i in idsStarSVMR])
idsGalaxySVMR = np.array([int(i) for i in idsGalaxySVMR])

#EC stands for ErrorCut
magGEC = -2.5*np.log10(matchedCat['cmodel_flux_g']/matchedCat['flux_zeromag_g'])
magREC = -2.5*np.log10(matchedCat['cmodel_flux_r']/matchedCat['flux_zeromag_r'])
magIEC = -2.5*np.log10(matchedCat['cmodel_flux_i']/matchedCat['flux_zeromag_i'])
magZEC = -2.5*np.log10(matchedCat['cmodel_flux_z']/matchedCat['flux_zeromag_z'])
magYEC = -2.5*np.log10(matchedCat['cmodel_flux_y']/matchedCat['flux_zeromag_y'])

ymagsEC = matchedCat['mag_y']
jmagsEC = matchedCat['mag_j']
hmagsEC = matchedCat['mag_h']
kmagsEC = matchedCat['mag_k']

chan1magsEC = matchedCat['mag_36']
chan2magsEC = matchedCat['mag_45']
chan3magsEC = matchedCat['mag_58']
chan4magsEC = matchedCat['mag_80']

def randomPlotColors(color1, color2, mag1Error, mag2Error, mag3Error, mag4Error, color1Deconv, color2Deconv, \
        color1Clustering, color2Clustering, xaxistitle, yaxistitle, filename):
    fig = plt.figure(figsize=(15, 15))

    ax1 = plt.subplot(321)

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

    plt.scatter(magSource1[galaxyIndicesLocal], magSource2[galaxyIndicesLocal], marker='o', c='blue', \
            label = 'Galaxies')
    plt.scatter(magSource1[starIndicesLocal], magSource2[starIndicesLocal], marker='o', c='red', \
            label='Stars')
    
    handles, labels = ax1.get_legend_handles_labels()
    # ax1.legend(handles, labels, loc='lower left')

    ax1.set_ylim([-2, 2])
    ax1.set_xlim([-2, 2])
    ax1.tick_params(axis='both', which='major', labelsize=12)
    ax1.tick_params(axis='both', which='minor', labelsize=12)
    ax1.set_xlabel(xaxistitle, fontdict={'fontsize': 12})
    ax1.set_ylabel(yaxistitle, fontdict={'fontsize': 12})
    #ax1.set_title('Distribution of Colors for All Sources')

    ax2 = plt.subplot(322)
    magSource1Bad = color1Deconv
    magSource2Bad = color2Deconv

    plt.scatter(magSource1Bad[deconvGalaxyIndices], magSource2Bad[deconvGalaxyIndices], c='blue', marker='o',\
            label = 'Galaxies')
    plt.scatter(magSource1Bad[deconvStarIndices], magSource2Bad[deconvStarIndices], c='red', marker='o',\
            label='Stars')

    ax2.set_ylim([-2, 2])
    ax2.set_xlim([-2, 2])
    ax2.tick_params(axis='both', which='major', labelsize=12)
    ax2.tick_params(axis='both', which='minor', labelsize=12)
    ax2.set_xlabel(xaxistitle, fontdict={'fontsize': 12})
    ax2.set_ylabel(yaxistitle, fontdict={'fontsize': 12})
    #ax2.set_title('Distribution of Colors for Misclassified Sources (XD)')
    handles, labels = ax2.get_legend_handles_labels()
    # ax2.legend(handles, labels, loc='lower left')

    # logistic regression IDs
    ax3 = plt.subplot(323)
    magSource1Bad = color1Clustering
    magSource2Bad = color2Clustering
    galaxyMask = ErrorCutIDsToMask(idsGalaxyLog)
    starMask = ErrorCutIDsToMask(idsStarLog)
    
    plt.scatter(magSource1Bad[galaxyMask], magSource2Bad[galaxyMask], c='blue', marker='o',\
            label = 'Galaxies')
    plt.scatter(magSource1Bad[starMask], magSource2Bad[starMask], c='red', marker='o',\
            label='Stars')
    ax3.set_ylim([-2, 2])
    ax3.set_xlim([-2, 2])
    ax3.tick_params(axis='both', which='major', labelsize=12)
    ax3.tick_params(axis='both', which='minor', labelsize=12)
    ax3.set_xlabel(xaxistitle, fontdict={'fontsize': 12})
    ax3.set_ylabel(yaxistitle, fontdict={'fontsize': 12})
    #ax3.set_title('Distribution of Colors for Misclassified Sources (LR)')
    handles, labels = ax3.get_legend_handles_labels()
    # ax3.legend(handles, labels, loc='lower left')

    # svm linear IDs
    ax4 = plt.subplot(324)
    magSource1Bad = color1Clustering
    magSource2Bad = color2Clustering
    galaxyMask = ErrorCutIDsToMask(idsGalaxySVML)
    starMask = ErrorCutIDsToMask(idsStarSVML)
    
    plt.scatter(magSource1Bad[galaxyMask], magSource2Bad[galaxyMask], c='blue', marker='o',\
            label = 'Galaxies')
    plt.scatter(magSource1Bad[starMask], magSource2Bad[starMask], c='red', marker='o',\
            label='Stars')
    ax4.set_ylim([-2, 2])
    ax4.set_xlim([-2, 2])
    ax4.tick_params(axis='both', which='major', labelsize=12)
    ax4.tick_params(axis='both', which='minor', labelsize=12)
    ax4.set_xlabel(xaxistitle, fontdict={'fontsize': 12})
    ax4.set_ylabel(yaxistitle, fontdict={'fontsize': 12})
    #ax4.set_title('Distribution of Colors for Misclassified Sources (SVM-Linear)')
    handles, labels = ax4.get_legend_handles_labels()
    # ax4.legend(handles, labels, loc='lower left')

    #SVM RBF IDs
    ax5 = plt.subplot(325)
    magSource1Bad = color1Clustering
    magSource2Bad = color2Clustering
    galaxyMask = ErrorCutIDsToMask(idsGalaxySVMR)
    starMask = ErrorCutIDsToMask(idsStarSVMR)
    
    plt.scatter(magSource1Bad[galaxyMask], magSource2Bad[galaxyMask], c='blue', marker='o',\
            label = 'Galaxies')
    plt.scatter(magSource1Bad[starMask], magSource2Bad[starMask], c='red', marker='o',\
            label='Stars')
    ax5.set_ylim([-2, 2])
    ax5.set_xlim([-2, 2])
    ax5.tick_params(axis='both', which='major', labelsize=12)
    ax5.tick_params(axis='both', which='minor', labelsize=12)
    ax5.set_xlabel(xaxistitle, fontdict={'fontsize': 12})
    ax5.set_ylabel(yaxistitle, fontdict={'fontsize': 12})
    #ax5.set_title('Distribution of Colors for Misclassified Sources (SVM-RBF)')
    handles, labels = ax5.get_legend_handles_labels()
    # ax5.legend(handles, labels, loc='lower left')


    # speed up performance
    fig.tight_layout()
    plt.savefig(filename)
    plt.close(fig)
    print 'done'



randomPlotColors(magGAll-magRAll, magRAll-magIAll, magGErrorAll, magRErrorAll, magRErrorAll, magIErrorAll, magGDeconv-magRDeconv, \
        magRDeconv-magIDeconv, magGEC-magREC, magREC-magIEC, 'g-r (AB)', 'r-i (AB)', 'data/SideBySide/gri.png')

randomPlotColors(magRAll-magIAll, magIAll-magZAll, magRErrorAll, magIErrorAll, magIErrorAll, magZErrorAll, \
        magRDeconv - magIDeconv, magIDeconv-magZDeconv, magREC - magIEC, magIEC - magZEC, \
        'r-i (AB)', 'i-z (AB)', 'data/SideBySide/riz.png')

randomPlotColors(magIAll-magZAll, magZAll-magYAll, magIErrorAll, magZErrorAll, magZErrorAll, magYErrorAll, \
        magIDeconv - magZDeconv, magZDeconv-magYDeconv, magIEC - magZEC, magZEC - magYEC, \
        'i-z (AB)', 'z-y (AB)', 'data/SideBySide/izy.png')

randomPlotColors(magZAll-magYAll, magYAll-jmagsAll, magZErrorAll, magYErrorAll, magYErrorAll, matchedCatAll['mag_j_error'], \
        magZDeconv - magYDeconv, magYDeconv-jmagsEC, magZEC - magYEC, magYEC - jmagsEC, \
        'z-y (AB)', 'y-J (AB)', 'data/SideBySide/zyj.png')

randomPlotColors(magYAll-jmagsAll, hmagsAll-kmagsAll, magYErrorAll, matchedCatAll['mag_j_error'], matchedCatAll['mag_h_error'] , \
        matchedCatAll['mag_k_error'], magYDeconv - jmagsDeconv, hmagsDeconv-kmagsDeconv, \
        magYEC - jmagsEC, hmagsEC - kmagsEC, \
        'y-J (AB)', 'H-K_s (AB)', 'data/SideBySide/yjhk.png')

randomPlotColors(jmagsAll-hmagsAll, hmagsAll-kmagsAll, matchedCatAll['mag_j_error'], matchedCatAll['mag_h_error'], \
        matchedCatAll['mag_h_error'], matchedCatAll['mag_k_error'], jmagsDeconv-hmagsDeconv, hmagsDeconv-kmagsDeconv, \
        jmagsEC-hmagsEC, hmagsEC-kmagsEC, 'J-H (AB)', 'H-K_s (AB)', 'data/SideBySide/jhk.png')

randomPlotColors(hmagsAll-kmagsAll, chan1magsAll-chan2magsAll, matchedCatAll['mag_h_error'], matchedCatAll['mag_k_error'], \
                matchedCatAll['mag_36error'], matchedCatAll['mag_45error'], hmagsDeconv - kmagsDeconv, chan1magsDeconv-chan2magsDeconv, \
                hmagsEC-kmagsEC, chan1magsEC-chan2magsEC, \
                r'$3.6\mu{}m - 4.5\mu{}m$', 'H-K_s', 'data/SideBySide/hkchan1chan2.png')
