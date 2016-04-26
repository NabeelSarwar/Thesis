import pyfits
import numpy as np
import extreme_deconvolution as xd
from sklearn.mixture import GMM
from sklearn.metrics import roc_curve, auc
from astroML.density_estimation import XDGMM

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

catalog = pyfits.getdata('MatchHellDeconv2.fits')

columns = ['magR', 'magI', 'magZ', 'magY',  'mag_j' ,'mag_h', 'mag_k', 'mag_36', 'mag_45', \
        'magRError', 'magIError', 'magZError', 'magYError', 'mag_j_error', \
        'mag_h_error', 'mag_k_error', 'mag_36error', 'mag_45error', 'mu_class', 'id']



cold = dict([(columns[i], i) for i in range(len(columns))])


catalog = np.array([catalog[columns[i]] for i in range(len(columns))]).transpose()

#goodIndices = []
#for i in range(len(catalog)):
   # if not np.any(np.isnan(catalog[i])):
    #    goodIndices.append(i)

#goodIndices = np.array(goodIndices)

#catalog = catalog[goodIndices]
starIndices = np.where(catalog[: , cold['mu_class']]==2)
galaxyIndices = np.where(catalog[: , cold['mu_class']]==1)

#the array index [0] is needed since np.where returns a tuple
PStar = 1.0 * len(starIndices[0]) / len(catalog)
PGalaxy = 1.0 * len(galaxyIndices[0])/ len(catalog)

print 'PStar is {0}'.format(PStar)

# can always get id
def makeEntry(index):
    colors = []
    data = catalog[index]
    colors.append(data[cold['magR']] - data[cold['magI']]) #0 
    colors.append(data[cold['magI']]- data[cold['magZ']]) #1
    colors.append(data[cold['magZ']] - data[cold['magY']]) #2 
    colors.append(data[cold['magY']] - data[cold['mag_j']]) #3
    colors.append(data[cold['mag_j']] - data[cold['mag_h']]) #4
    colors.append(data[cold['mag_h']] - data[cold['mag_k']]) #5
    colors.append(data[cold['mag_k']] - data[cold['mag_36']]) #6
    colors.append(data[cold['mag_36']] - data[cold['mag_45']]) #7
    colors = np.array(colors)
    if np.any(np.isnan(colors)):
        print index
    return colors 

def  makeNoiseMatrix(index):
    noise = np.zeros(64).reshape(8, 8)

    data = catalog[index]

    noise[0, 0] = data[cold['magRError']]**2 + data[cold['magIError']]**2
    noise[1, 1] = data[cold['magIError']]**2 + data[cold['magZError']]**2
    noise[2, 2] = data[cold['magZError']]**2 + data[cold['magYError']]**2
    noise[3, 3] = data[cold['magYError']]**2 + data[cold['mag_j_error']]**2
    noise[4, 4] = data[cold['mag_j_error']]**2 + data[cold['mag_h_error']]**2
    noise[5, 5] = data[cold['mag_h_error']]**2 + data[cold['mag_k_error']]**2
    noise[6, 6] = data[cold['mag_k_error']]**2 + data[cold['mag_36error']]**2
    noise[7, 7] = data[cold['mag_36error']]**2 + data[cold['mag_45error']]**2

    noise[0, 1] = -data[cold['magIError']]**2
    noise[1, 2] = -data[cold['magZError']]**2
    noise[2, 3] = -data[cold['magYError']]**2
    noise[3, 4] = -data[cold['mag_j_error']]**2
    noise[4, 5] = -data[cold['mag_h_error']]**2
    noise[5, 6] = -data[cold['mag_k_error']]**2
    noise[6, 7] = -data[cold['mag_36error']]**2

    noise[1, 0] = noise[0, 1]
    noise[2, 1] = noise[1, 2]
    noise[3, 2] = noise[2, 3]
    noise[4, 3] = noise[3, 4]
    noise[5, 4] = noise[4, 5]
    noise[6, 5] = noise[5, 6]
    noise[7, 6] = noise[6, 7]

    return noise



# 1 is star prediction, 0 is no star
def predictStarDouble(X, Xerr, index):
    if np.any(np.isnan(Xerr)):
        print index
    #numerator = PStar * np.exp(clfstar.logprob_a(X, Xerr))
    #demominator = PStar * np.exp(clfstar.logprob_a(X, Xerr)) + PGalaxy * np.exp(clfgalaxy.logprob_a(X, Xerr))
    #P(Star|X, XErr)

    fraction = np.log(PStar) + clfstar.logprob_a(X, Xerr)[0] \
            - np.logaddexp(np.log(PStar) +  clfstar.logprob_a(X, Xerr)[0], np.log(PGalaxy) + clfgalaxy.logprob_a(X, Xerr)[0])
    
    fraction = np.sum(np.exp(fraction))

    if np.isnan(fraction):
        raise Exception('Invalid Fractions')

    return fraction

def predictStar(X, Xerr, index):
    fraction = predictStarDouble(X, Xerr, index)
    if fraction >= 0.5:
        return 1
    else:
        return 0

def generateReport(predictions, results):
    truepositive = 0
    falsepositive = 0
    truenegative = 0
    falsenegative = 0

    if len(results) != len(predictions):
        print 'Incorrect predictions'
    for i in range(len(predictions)):
        pred = predictions[i]
        actual = results[i]
        if actual == 1:
            if pred == 1:
                truepositive +=1
            else:
                falsenegative +=1
        elif actual==0:
            if pred==1:
                falsepositive+=1
            else:
                truenegative+=1
        else:
            print 'Bad result'
            print actual
    report = {}
    report['TPR'] = 1.0*truepositive/len(results)
    report['TNR'] = 1.0*truenegative/len(results)
    report['FPR'] = 1.0*falsepositive/len(results)
    report['FNR'] = 1.0*falsenegative/len(results)
    report['Precision'] = 1.0 * truepositive / (truepositive + falsepositive)
    report['Recall'] = 1.0 * truepositive/(truepositive+falsenegative)
    report['Accuracy'] = 1.0 * (truepositive + truenegative) / len(results)
    return report

def writeReport(report, fileName):
    buf = open(fileName, 'w')
    buf.write('TPR: {0}\n'.format(report['TPR']))
    buf.write('TNR: {0}\n'.format(report['TNR']))
    buf.write('FPR: {0}\n'.format(report['FPR']))
    buf.write('FNR: {0}\n'.format(report['FNR']))
    buf.write('Precision: {0}\n'.format(report['Precision']))
    buf.write('Recall: {0}\n'.format(report['Recall']))
    buf.write('Accuracy: {0}\n'.format(report['Accuracy']))
    buf.close()

def generateROCCurve(predictions, results):
    fig = plt.figure()
    fpr, tpr, thresholds = roc_curve(results, predictions)
    plt.plot(fpr, tpr, label='AUC={0}'.format(auc(fpr, tpr)))
    plt.plot([0, 1], [0,1], 'k--')
    plt.xlim([0, 1])
    plt.ylim([0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Curve of Point Source Classification')
    plt.legend(loc="lower right")
    plt.savefig('data/PStar/ROC.png')

starTrainNumber = int(np.floor(len(starIndices[0]) * 0.8))
starTestNumber = len(starIndices[0]) - starTrainNumber
print 'Number of test stars: {0}'.format(starTestNumber)
starIndicesTrain = np.random.choice(starIndices[0], starTrainNumber, replace=False)
starIndicesTest = []

for index in starIndices[0]:
    if not index in starIndicesTrain:
        starIndicesTest.append(index)
starIndicesTest = np.array(starIndicesTest)

if len(starIndicesTest) != starTestNumber:
    raise Exception('False Star Test numbers')

galaxyTrainNumber = int(np.floor(len(galaxyIndices[0]) * 0.8))
galaxyTestNumber = len(galaxyIndices[0])- galaxyTrainNumber

print 'Number of test galaxies: {0}'.format(galaxyTestNumber)

galaxyIndicesTrain = np.random.choice(galaxyIndices[0], galaxyTrainNumber, replace=False)
galaxyIndicesTest = []
for index in galaxyIndices[0]:
    if not index in galaxyIndicesTrain:
        galaxyIndicesTest.append(index)
galaxyIndicesTest = np.array(galaxyIndicesTest)

if len(galaxyIndicesTest) != galaxyTestNumber:
    raise Exception('False Galaxy Test numbers')

nGaussiansStar = 19
nGaussiansGalaxy = 19

print 'Making Arrays'
print 'Making X'
XTrainStar = np.array([makeEntry(index) for index in starIndicesTrain])
XTestStar = np.array([makeEntry(index) for index in starIndicesTest])
XTrainGalaxy = np.array([makeEntry(index) for index in galaxyIndicesTrain])
XTestGalaxy = np.array([makeEntry(index) for index in galaxyIndicesTest])

print 'Making Errors'
XErrTrainStar = np.array([makeNoiseMatrix(index) for index in starIndicesTrain])
XErrTestStar = np.array([makeNoiseMatrix(index) for index in starIndicesTest])
XErrTrainGalaxy = np.array([makeNoiseMatrix(index) for index in galaxyIndicesTrain])
XErrTestGalaxy = np.array([makeNoiseMatrix(index) for index in galaxyIndicesTest])

#convolving
print 'Estimating Gaussians'
GMMStar = GMM(nGaussiansStar, n_iter = 10, covariance_type='full').fit(XTrainStar)
GMMGalaxy = GMM(nGaussiansGalaxy, n_iter=10, covariance_type='full').fit(XTrainGalaxy)

ampstar = GMMStar.weights_
meanstar = GMMStar.means_
covarstar = GMMStar.covars_

ampgalaxy = GMMGalaxy.weights_
meangalaxy = GMMGalaxy.means_
covargalaxy = GMMGalaxy.covars_


# Results are saved in `amp`, `mean`, and `covar`
print 'Deconvolving star'

xd.extreme_deconvolution(XTrainStar, XErrTrainStar, ampstar, meanstar, covarstar)

clfstar = XDGMM(nGaussiansStar)
clfstar.alpha = ampstar
clfstar.mu = meanstar
clfstar.V = covarstar


print 'Deconvolving galaxies'
xd.extreme_deconvolution(XTrainGalaxy, XErrTrainGalaxy, ampgalaxy, meangalaxy, covargalaxy)

clfgalaxy = XDGMM(nGaussiansGalaxy)
clfgalaxy.alpha = ampgalaxy
clfgalaxy.mu = meangalaxy
clfgalaxy.V = covargalaxy

print 'Predicting'
# need to pass XTestStar[i] and XTestGalaxy[i] as np.array([XTestStar[i]]) because internally it assumes 2D matrix
starPredictions = np.array([predictStar(np.array([XTestStar[i]]), np.array([XErrTestStar[i]]), starIndicesTest[i]) for i in range(starTestNumber)])

galaxyPredictions = np.array([predictStar(np.array([XTestGalaxy[i]]), np.array([XErrTestGalaxy[i]]), galaxyIndicesTest[i]) for i in range(galaxyTestNumber)])

predictions = np.array(starPredictions.tolist() +  galaxyPredictions.tolist())
results = np.array([1 for i in range(len(starPredictions))] + [0 for i in range(len(galaxyPredictions))])

report = generateReport(predictions, results)
writeReport(report, 'data/deconv/deconvresults.txt')

generateROCCurve(predictions, results)

#analyze bad classifications
badStarPredictions = np.where(starPredictions == 0)[0]
badStarIndices = starIndicesTest[badStarPredictions]

# these are indices of the array where there are bad predictions
# use these to get the bad indices in the catalog
badGalaxyPredictions = np.where(galaxyPredictions == 1)[0]
badGalaxyIndices = galaxyIndicesTest[badGalaxyPredictions]

starPredictionsNumbers = np.array([predictStarDouble(np.array([XTestStar[i]]), np.array([XErrTestStar[i]]), starIndicesTest[i]) for i in range(starTestNumber)])

galaxyPredictionsNumbers = np.array([predictStarDouble(np.array([XTestGalaxy[i]]), np.array([XErrTestGalaxy[i]]), galaxyIndicesTest[i]) for i in range(galaxyTestNumber)])

# these are the bad probabilities
badStarPredictionsDeconv = starPredictionsNumbers[badStarPredictions]
badGalaxyPredictionsDeconv = galaxyPredictionsNumbers[badGalaxyPredictions]

starSaveInformation = np.array([badStarIndices, badStarPredictions]).transpose()
galaxySaveInformation = np.array([badGalaxyIndices, badGalaxyPredictions]).transpose()

np.savetxt('data/deconv/starBadPredictionIndexProbability.txt', starSaveInformation)
np.savetxt('data/deconv/galaxyBadPredictionIndexProbability.txt', galaxySaveInformation)


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
chan3magsEC = matchedCatEC['mag_58']
chan4magsEC = matchedCatEC['mag_80']

# indices correspond to indices in catalog above
def AnalyzeBadPoints(predictions, indices, magnitudes, xaxistitle):
    color = magnitudes[indices]
    plt.ylabel(r'$P(Star)$')
    plt.xlabel(xaxistitle)

    plt.scatter(color, predictions)
    # do it for x

    xticks, xticklabels = plt.xticks()

    xmin = (3*xticks[0] - xticks[1])/2.
    xmax = (3*xticks[-1] - xticks[-2])/2.

    plt.xlim(xmin, xmax)
    plt.xticks(xticks) 

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

        ax2 = plt.subplot(222)
        indices2 = ErrorCutIDsToMask(matchedCatEC, indices2)
        AnalyzeBadPoints(predictions2, indices2, magnitudes2, xaxistitle)

        ax3 = plt.subplot(223)
        indices3 = ErrorCutIDsToMask(matchedCatEC, indices3)
        AnalyzeBadPoints(predictions3, indices3, magnitudes2, xaxistitle)
        
        ax4 = plt.subplot(224)
        indices4 = ErrorCutIDsToMask(matchedCatEC, indices4)
        AnalyzeBadPoints(predictions4, indices4, magnitudes2, xaxistitle)
        plt.savefig(filename)


# DO THE SIDE BY SIDE ALL PSTARS
MultipleViewAnalyzeBadPoints(matchedCatEC, badStarPredictions, badStarIndices, catalog[:, cold['magR']], \
        magREC, starLogPredictions, idsStarLog, \
        starSVMLPredictions, idsStarSVML, \
        starSVMRPredictions, idsStarSVMR, 
        'r (AB)', 'data/PStar/side_false_negatives_r_star.png')
MultipleViewAnalyzeBadPoints(matchedCatEC, badStarPredictions, badStarIndices, catalog[:, cold['magI']], \
        magIEC, starLogPredictions, idsStarLog, \
        starSVMLPredictions, idsStarSVML, \
        starSVMRPredictions, idsStarSVMR, 
        'i (AB)', 'data/PStar/side_false_negatives_i_star.png')
MultipleViewAnalyzeBadPoints(matchedCatEC, badStarPredictions, badStarIndices, catalog[:, cold['magZ']], \
        magZEC, starLogPredictions, idsStarLog, \
        starSVMLPredictions, idsStarSVML, \
        starSVMRPredictions, idsStarSVMR, 
        'z (AB)', 'data/PStar/side_false_negatives_z_star.png')

MultipleViewAnalyzeBadPoints(matchedCatEC, badStarPredictions, badStarIndices, catalog[:, cold['magY']], \
        magYEC, starLogPredictions, idsStarLog, \
        starSVMLPredictions, idsStarSVML, \
        starSVMRPredictions, idsStarSVMR, 
        'y (AB)', 'data/PStar/side_false_negatives_y_star.png')

MultipleViewAnalyzeBadPoints(matchedCatEC, badStarPredictions, badStarIndices, catalog[:, cold['mag_j']], \
        jmagsEC, starLogPredictions, idsStarLog, \
        starSVMLPredictions, idsStarSVML, \
        starSVMRPredictions, idsStarSVMR, 
        'J (AB)', 'data/PStar/side_false_negatives_j_star.png')

MultipleViewAnalyzeBadPoints(matchedCatEC, badStarPredictions, badStarIndices, catalog[:, cold['mag_h']], \
        hmagsEC, starLogPredictions, idsStarLog, \
        starSVMLPredictions, idsStarSVML, \
        starSVMRPredictions, idsStarSVMR, 
        'H (AB)', 'data/PStar/side_false_negatives_h_star.png')

MultipleViewAnalyzeBadPoints(matchedCatEC, badStarPredictions, badStarIndices, catalog[:, cold['mag_k']], \
        kmagsEC, starLogPredictions, idsStarLog, \
        starSVMLPredictions, idsStarSVML, \
        starSVMRPredictions, idsStarSVMR, 
        'K_s (AB)', 'data/PStar/side_false_negatives_k_star.png')

MultipleViewAnalyzeBadPoints(matchedCatEC, badStarPredictions, badStarIndices, catalog[:, cold['mag_36']], \
        chan1magsEC, starLogPredictions, idsStarLog, \
        starSVMLPredictions, idsStarSVML, \
        starSVMRPredictions, idsStarSVMR, 
        r'$3.6\mu{}m$ (AB)', 'data/PStar/side_false_negatives_chan1_star.png')
MultipleViewAnalyzeBadPoints(matchedCatEC, badStarPredictions, badStarIndices, catalog[:, cold['mag_45']], \
        chan2magsEC, starLogPredictions, idsStarLog, \
        starSVMLPredictions, idsStarSVML, \
        starSVMRPredictions, idsStarSVMR, 
        r'$4.5\mu{}m$ (AB)', 'data/PStar/side_false_negatives_chan2_star.png')

MultipleViewAnalyzeBadPoints(matchedCatEC, badGalaxyPredictions, badGalaxyIndices, catalog[:, cold['magR']], \
        magREC, galaxyLogPredictions, idsGalaxyLog, \
        galaxySVMLPredictions, idsGalaxySVML, \
        galaxySVMRPredictions, idsGalaxySVMR, 
        'r (AB)', 'data/PStar/side_false_positives_r_galaxy.png')
MultipleViewAnalyzeBadPoints(matchedCatEC, badGalaxyPredictions, badGalaxyIndices, catalog[:, cold['magI']], \
        magIEC, galaxyLogPredictions, idsGalaxyLog, \
        galaxySVMLPredictions, idsGalaxySVML, \
        galaxySVMRPredictions, idsGalaxySVMR, 
        'i (AB)', 'data/PStar/side_false_positives_i_galaxy.png')
MultipleViewAnalyzeBadPoints(matchedCatEC, badGalaxyPredictions, badGalaxyIndices, catalog[:, cold['magZ']], \
        magZEC, galaxyLogPredictions, idsGalaxyLog, \
        galaxySVMLPredictions, idsGalaxySVML, \
        galaxySVMRPredictions, idsGalaxySVMR, 
        'z (AB)', 'data/PStar/side_false_positives_z_galaxy.png')

MultipleViewAnalyzeBadPoints(matchedCatEC, badGalaxyPredictions, badGalaxyIndices, catalog[:, cold['magY']], \
        magYEC, galaxyLogPredictions, idsGalaxyLog, \
        galaxySVMLPredictions, idsGalaxySVML, \
        galaxySVMRPredictions, idsGalaxySVMR, 
        'y (AB)', 'data/PStar/side_false_positives_y_galaxy.png')

MultipleViewAnalyzeBadPoints(matchedCatEC, badGalaxyPredictions, badGalaxyIndices, catalog[:, cold['mag_j']], \
        jmagsEC, galaxyLogPredictions, idsGalaxyLog, \
        galaxySVMLPredictions, idsGalaxySVML, \
        galaxySVMRPredictions, idsGalaxySVMR, 
        'J (AB)', 'data/PStar/side_false_positives_j_galaxy.png')

MultipleViewAnalyzeBadPoints(matchedCatEC, badGalaxyPredictions, badGalaxyIndices, catalog[:, cold['mag_h']], \
        hmagsEC, galaxyLogPredictions, idsGalaxyLog, \
        galaxySVMLPredictions, idsGalaxySVML, \
        galaxySVMRPredictions, idsGalaxySVMR, 
        'H (AB)', 'data/PStar/side_false_positives_h_galaxy.png')

MultipleViewAnalyzeBadPoints(matchedCatEC, badGalaxyPredictions, badGalaxyIndices, catalog[:, cold['mag_k']], \
        kmagsEC, galaxyLogPredictions, idsGalaxyLog, \
        galaxySVMLPredictions, idsGalaxySVML, \
        galaxySVMRPredictions, idsGalaxySVMR, 
        'K_s (AB)', 'data/PStar/side_false_positives_k_galaxy.png')

MultipleViewAnalyzeBadPoints(matchedCatEC, badGalaxyPredictions, badGalaxyIndices, catalog[:, cold['mag_36']], \
        chan1magsEC, galaxyLogPredictions, idsGalaxyLog, \
        galaxySVMLPredictions, idsGalaxySVML, \
        galaxySVMRPredictions, idsGalaxySVMR, 
        r'$3.6\mu{}m$ (AB)', 'data/PStar/side_false_positives_chan1_galaxy.png')
MultipleViewAnalyzeBadPoints(matchedCatEC, badGalaxyPredictions, badGalaxyIndices, catalog[:, cold['mag_45']], \
        chan2magsEC, galaxyLogPredictions, idsGalaxyLog, \
        galaxySVMLPredictions, idsGalaxySVML, \
        galaxySVMRPredictions, idsGalaxySVMR, 
        r'$4.5\mu{}m$ (AB)', 'data/PStar/side_false_positives_chan2_galaxy.png')
