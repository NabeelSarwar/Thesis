import pyfits
import numpy as np
import extreme_deconvolution as xd
from sklearn.mixture import GMM
from sklearn.metrics import roc_curve, auc
from astroML.density_estimation import XDGMM
import scipy.misc as misc

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

catalog = pyfits.getdata('MatchHellDeconv2.fits')

columns = ['magR', 'magI', 'magZ', 'magY',  'mag_j' ,'mag_h', 'mag_k', 'mag_36', 'mag_45', \
        'magRError', 'magIError', 'magZError', 'magYError', 'mag_j_error', \
        'mag_h_error', 'mag_k_error', 'mag_36error', 'mag_45error', 'mu_class', 'id']



cold = dict([(columns[i], i) for i in range(len(columns))])


catalog = np.array([catalog[columns[i]] for i in range(len(columns))]).transpose()

badStarIndices = np.genfromtxt('data/deconv/starBadPredictionIndexProbability.txt')[:, 0]
badGalaxyIndices = np.genfromtxt('data/deconv/galaxyBadPredictionIndexProbability.txt')[:, 0]
badStarIndices = np.array([int(i) for i in badStarIndices])
badGalaxyIndices = np.array([int(i) for i in badGalaxyIndices])
badStarPredictionsDeconv  = np.genfromtxt('data/deconv/starBadPredictionIndexProbability.txt')[:, 1]
badGalaxyPredictionsDeconv = np.genfromtxt('data/deconv/galaxyBadPredictionIndexProbability.txt')[:, 1]

matchedCatEC = pyfits.getdata('MatchHellErrorCut2.fits')
idsStarLog = np.genfromtxt('data/logistic/badStarIDs.txt')[:, 0]
idsGalaxyLog = np.genfromtxt('data/logistic/badGalaxyIDs.txt')[:, 0]
idsStarLog = np.array([int(i) for i in idsStarLog])
idsGalaxyLog = np.array([int(i) for i in idsGalaxyLog])
starLogPredictions  = np.genfromtxt('data/logistic/badStarIDs.txt')[:, 1]
galaxyLogPredictions = np.genfromtxt('data/logistic/badGalaxyIDs.txt')[:, 1]

idsStarSVML = np.genfromtxt('data/svm/badStarIDslinear.txt')[:, 0]
idsGalaxySVML = np.genfromtxt('data/svm/badGalaxyIDslinear.txt')[:, 0]
idsStarSVML = np.array([int(i) for i in idsStarSVML])
idsGalaxySVML = np.array([int(i) for i in idsGalaxySVML])
starSVMLPredictions = np.genfromtxt('data/svm/badStarIDslinear.txt')[:, 1]
galaxySVMLPredictions = np.genfromtxt('data/svm/badGalaxyIDslinear.txt')[:, 1]

idsStarSVMR = np.genfromtxt('data/svm/badStarIDsRBF.txt')[:, 0]
idsGalaxySVMR = np.genfromtxt('data/svm/badGalaxyIDsRBF.txt')[:, 0]
idsStarSVMR = np.array([int(i) for i in idsStarSVMR])
idsGalaxySVMR = np.array([int(i) for i in idsGalaxySVMR])
starSVMRPredictions = np.genfromtxt('data/svm/badStarIDsRBF.txt')[:, 1]
galaxySVMRPredictions = np.genfromtxt('data/svm/badGalaxyIDsRBF.txt')[:, 1]

#EC stands for ErrorCut
magGEC = -2.5*np.log10(matchedCatEC['cmodel_flux_g']/matchedCatEC['flux_zeromag_g'])
magREC = -2.5*np.log10(matchedCatEC['cmodel_flux_r']/matchedCatEC['flux_zeromag_r'])
magIEC = -2.5*np.log10(matchedCatEC['cmodel_flux_i']/matchedCatEC['flux_zeromag_i'])
magZEC = -2.5*np.log10(matchedCatEC['cmodel_flux_z']/matchedCatEC['flux_zeromag_z'])
magYEC = -2.5*np.log10(matchedCatEC['cmodel_flux_y']/matchedCatEC['flux_zeromag_y'])

ymagsEC = matchedCatEC['mag_y']
jmagsEC = matchedCatEC['mag_j']
hmagsEC = matchedCatEC['mag_h']
kmagsEC = matchedCatEC['mag_k']

chan1magsEC = matchedCatEC['mag_36']
chan2magsEC = matchedCatEC['mag_45']
chan4magsEC = matchedCatEC['mag_58']
chan4magsEC = matchedCatEC['mag_80']

# indices correspond to indices in catalog above
def AnalyzeBadPoints(predictions, indices, magnitudes, xaxistitle):
    color = magnitudes[indices]
    plt.ylabel(r'$P(Star)$')
    plt.xlabel(xaxistitle)

    plt.scatter(color, predictions)

    
# to analyze bad SVM points
def ErrorCutIDsToMask(matchedCat, idArray):
    indices = []
    column = matchedCat['id']
    for i in range(len(column)):
        if column[i] in idArray:
            indices.append(i)
    mask = np.array([True if i in indices else False for i in range(len(column))])
    return mask


# 3 of the methods share the same magnitudes
def MultipleViewAnalyzeBadPoints(matchedCatEC, predictions1, indices1, magnitudes1, \
        magnitudes2, predictions2, indices2, \
        predictions3, indices3, \
        predictions4, indices4, \
        xaxistitle, filename):

        plt.figure(figsize=(15, 15))
        print xaxistitle

        ax1 = plt.subplot(221)
        AnalyzeBadPoints(predictions1, indices1, magnitudes1, xaxistitle)
        ax1.set_xlim(18, 28)

        ax2 = plt.subplot(222)
        indices2 = ErrorCutIDsToMask(matchedCatEC, indices2)
        AnalyzeBadPoints(predictions2, indices2, magnitudes2, xaxistitle)
        ax2.set_xlim(18, 28)

        ax3 = plt.subplot(223)
        indices3 = ErrorCutIDsToMask(matchedCatEC, indices3)
        AnalyzeBadPoints(predictions3, indices3, magnitudes2, xaxistitle)
        ax3.set_xlim(18, 28)
       
        ax4 = plt.subplot(224)
        indices4 = ErrorCutIDsToMask(matchedCatEC, indices4)
        AnalyzeBadPoints(predictions4, indices4, magnitudes2, xaxistitle)
        ax4.set_xlim(18, 28)

        plt.savefig(filename)

# DO THE SIDE BY SIDE ALL PSTARS
MultipleViewAnalyzeBadPoints(matchedCatEC, badStarPredictionsDeconv, badStarIndices, catalog[:, cold['magR']], \
        magREC, starLogPredictions, idsStarLog, \
        starSVMLPredictions, idsStarSVML, \
        starSVMRPredictions, idsStarSVMR, 
        'r (AB)', 'data/PStar/side_false_negatives_r_star.png')
MultipleViewAnalyzeBadPoints(matchedCatEC, badStarPredictionsDeconv, badStarIndices, catalog[:, cold['magI']], \
        magIEC, starLogPredictions, idsStarLog, \
        starSVMLPredictions, idsStarSVML, \
        starSVMRPredictions, idsStarSVMR, 
        'i (AB)', 'data/PStar/side_false_negatives_i_star.png')
MultipleViewAnalyzeBadPoints(matchedCatEC, badStarPredictionsDeconv, badStarIndices, catalog[:, cold['magZ']], \
        magZEC, starLogPredictions, idsStarLog, \
        starSVMLPredictions, idsStarSVML, \
        starSVMRPredictions, idsStarSVMR, 
        'z (AB)', 'data/PStar/side_false_negatives_z_star.png')

MultipleViewAnalyzeBadPoints(matchedCatEC, badStarPredictionsDeconv, badStarIndices, catalog[:, cold['magY']], \
        magYEC, starLogPredictions, idsStarLog, \
        starSVMLPredictions, idsStarSVML, \
        starSVMRPredictions, idsStarSVMR, 
        'y (AB)', 'data/PStar/side_false_negatives_y_star.png')

MultipleViewAnalyzeBadPoints(matchedCatEC, badStarPredictionsDeconv, badStarIndices, catalog[:, cold['mag_j']], \
        jmagsEC, starLogPredictions, idsStarLog, \
        starSVMLPredictions, idsStarSVML, \
        starSVMRPredictions, idsStarSVMR, 
        'J (AB)', 'data/PStar/side_false_negatives_j_star.png')

MultipleViewAnalyzeBadPoints(matchedCatEC, badStarPredictionsDeconv, badStarIndices, catalog[:, cold['mag_h']], \
        hmagsEC, starLogPredictions, idsStarLog, \
        starSVMLPredictions, idsStarSVML, \
        starSVMRPredictions, idsStarSVMR, 
        'H (AB)', 'data/PStar/side_false_negatives_h_star.png')

MultipleViewAnalyzeBadPoints(matchedCatEC, badStarPredictionsDeconv, badStarIndices, catalog[:, cold['mag_k']], \
        kmagsEC, starLogPredictions, idsStarLog, \
        starSVMLPredictions, idsStarSVML, \
        starSVMRPredictions, idsStarSVMR, 
        'K_s (AB)', 'data/PStar/side_false_negatives_k_star.png')

MultipleViewAnalyzeBadPoints(matchedCatEC, badStarPredictionsDeconv, badStarIndices, catalog[:, cold['mag_36']], \
        chan1magsEC, starLogPredictions, idsStarLog, \
        starSVMLPredictions, idsStarSVML, \
        starSVMRPredictions, idsStarSVMR, 
        r'$3.6\mu{}m$ (AB)', 'data/PStar/side_false_negatives_chan1_star.png')
MultipleViewAnalyzeBadPoints(matchedCatEC, badStarPredictionsDeconv, badStarIndices, catalog[:, cold['mag_45']], \
        chan2magsEC, starLogPredictions, idsStarLog, \
        starSVMLPredictions, idsStarSVML, \
        starSVMRPredictions, idsStarSVMR, 
        r'$4.5\mu{}m$ (AB)', 'data/PStar/side_false_negatives_chan2_star.png')

MultipleViewAnalyzeBadPoints(matchedCatEC, badGalaxyPredictionsDeconv, badGalaxyIndices, catalog[:, cold['magR']], \
        magREC, galaxyLogPredictions, idsGalaxyLog, \
        galaxySVMLPredictions, idsGalaxySVML, \
        galaxySVMRPredictions, idsGalaxySVMR, 
        'r (AB)', 'data/PStar/side_false_positives_r_galaxy.png')
MultipleViewAnalyzeBadPoints(matchedCatEC, badGalaxyPredictionsDeconv, badGalaxyIndices, catalog[:, cold['magI']], \
        magIEC, galaxyLogPredictions, idsGalaxyLog, \
        galaxySVMLPredictions, idsGalaxySVML, \
        galaxySVMRPredictions, idsGalaxySVMR, 
        'i (AB)', 'data/PStar/side_false_positives_i_galaxy.png')
MultipleViewAnalyzeBadPoints(matchedCatEC, badGalaxyPredictionsDeconv, badGalaxyIndices, catalog[:, cold['magZ']], \
        magZEC, galaxyLogPredictions, idsGalaxyLog, \
        galaxySVMLPredictions, idsGalaxySVML, \
        galaxySVMRPredictions, idsGalaxySVMR, 
        'z (AB)', 'data/PStar/side_false_positives_z_galaxy.png')

MultipleViewAnalyzeBadPoints(matchedCatEC, badGalaxyPredictionsDeconv, badGalaxyIndices, catalog[:, cold['magY']], \
        magYEC, galaxyLogPredictions, idsGalaxyLog, \
        galaxySVMLPredictions, idsGalaxySVML, \
        galaxySVMRPredictions, idsGalaxySVMR, 
        'y (AB)', 'data/PStar/side_false_positives_y_galaxy.png')

MultipleViewAnalyzeBadPoints(matchedCatEC, badGalaxyPredictionsDeconv, badGalaxyIndices, catalog[:, cold['mag_j']], \
        jmagsEC, galaxyLogPredictions, idsGalaxyLog, \
        galaxySVMLPredictions, idsGalaxySVML, \
        galaxySVMRPredictions, idsGalaxySVMR, 
        'J (AB)', 'data/PStar/side_false_positives_j_galaxy.png')

MultipleViewAnalyzeBadPoints(matchedCatEC, badGalaxyPredictionsDeconv, badGalaxyIndices, catalog[:, cold['mag_h']], \
        hmagsEC, galaxyLogPredictions, idsGalaxyLog, \
        galaxySVMLPredictions, idsGalaxySVML, \
        galaxySVMRPredictions, idsGalaxySVMR, 
        'H (AB)', 'data/PStar/side_false_positives_h_galaxy.png')

MultipleViewAnalyzeBadPoints(matchedCatEC, badGalaxyPredictionsDeconv, badGalaxyIndices, catalog[:, cold['mag_k']], \
        kmagsEC, galaxyLogPredictions, idsGalaxyLog, \
        galaxySVMLPredictions, idsGalaxySVML, \
        galaxySVMRPredictions, idsGalaxySVMR, 
        'K_s (AB)', 'data/PStar/side_false_positives_k_galaxy.png')

MultipleViewAnalyzeBadPoints(matchedCatEC, badGalaxyPredictionsDeconv, badGalaxyIndices, catalog[:, cold['mag_36']], \
        chan1magsEC, galaxyLogPredictions, idsGalaxyLog, \
        galaxySVMLPredictions, idsGalaxySVML, \
        galaxySVMRPredictions, idsGalaxySVMR, 
        r'$3.6\mu{}m$ (AB)', 'data/PStar/side_false_positives_chan1_galaxy.png')
MultipleViewAnalyzeBadPoints(matchedCatEC, badGalaxyPredictionsDeconv, badGalaxyIndices, catalog[:, cold['mag_45']], \
        chan2magsEC, galaxyLogPredictions, idsGalaxyLog, \
        galaxySVMLPredictions, idsGalaxySVML, \
        galaxySVMRPredictions, idsGalaxySVMR, 
        r'$4.5\mu{}m$ (AB)', 'data/PStar/side_false_positives_chan2_galaxy.png')
