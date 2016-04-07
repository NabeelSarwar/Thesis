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

# can always get id
def makeEntry(index):
    colors = []
    data = catalog[index]
    colors.append(data[cold['magR']] - data[cold['magI']])
    colors.append(data[cold['magI']]- data[cold['magZ']])
    colors.append(data[cold['magZ']] - data[cold['magY']])
    colors.append(data[cold['mag_j']] - data[cold['mag_h']])
    colors.append(data[cold['mag_h']] - data[cold['mag_k']])
    colors.append(data[cold['mag_k']] - data[cold['mag_36']])
    colors.append(data[cold['mag_36']] - data[cold['mag_45']])
    colors = np.array(colors)

    if np.any(np.isnan(colors)):
        print index
    return colors 

def  makeNoiseMatrix(index):
    noise = np.zeros(49).reshape(7, 7)

    data = catalog[index]

    noise[0, 0] = data[cold['magRError']]**2 + data[cold['magIError']]**2
    noise[1, 1] = data[cold['magIError']]**2 + data[cold['magZError']]**2
    noise[2, 2] = data[cold['magZError']]**2 + data[cold['magYError']]**2

    noise[3, 3] = data[cold['mag_j_error']]**2 + data[cold['mag_h_error']]**2
    noise[4, 4] = data[cold['mag_h_error']]**2 + data[cold['mag_k_error']]**2
    noise[5, 5] = data[cold['mag_k_error']]**2 + data[cold['mag_36error']]**2
    noise[6, 6] = data[cold['mag_36error']]**2 + data[cold['mag_45error']]**2

    noise[0, 1] = -data[cold['magIError']]**2
    noise[1, 2] = -data[cold['magZError']]**2
    noise[3, 4] = -data[cold['mag_h_error']]**2
    noise[4, 5] = -data[cold['mag_k_error']]**2
    noise[5, 6] = -data[cold['mag_36error']]**2

    noise[1, 0] = noise[0, 1]
    noise[2, 1] = noise[1, 2]
    noise[4, 3] = noise[3, 4]
    noise[5, 4] = noise[4, 5]
    noise[6, 5] = noise[5, 6]

    return noise


# 1 is star prediction, 0 is no star
def predictStarDouble(X, Xerr, index):
    if np.any(np.isnan(Xerr)):
        print index
    #numerator = PStar * np.exp(clfstar.logL(X, Xerr))
    #demominator = PStar * np.exp(clfstar.logL(X, Xerr)) + PGalaxy * np.exp(clfgalaxy.logL(X, Xerr))
    #P(Star|X, XErr)

    fraction = np.log(PStar) + clfstar.logL(X, Xerr) \
            - np.logaddexp(np.log(PStar) +  clfstar.logL(X, Xerr), np.log(PGalaxy) + clfgalaxy.logL(X, Xerr))

    if np.isnan(fraction):
        raise Exception('Invalid Fractions')

    fraction = np.exp(fraction)
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
galaxyIndicesTrain = np.random.choice(galaxyIndices[0], galaxyTrainNumber, replace=False)
galaxyIndicesTest = []
for index in galaxyIndices[0]:
    if not index in galaxyIndicesTrain:
        galaxyIndicesTest.append(index)
galaxyIndicesTest = np.array(galaxyIndicesTest)

if len(galaxyIndicesTest) != galaxyTestNumber:
    raise Exception('False Galaxy Test numbers')

nGaussiansStar = 13
nGaussiansGalaxy = 27

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
badStarPredictions = starPredictionsNumbers[badStarPredictions]
badGalaxyPredictions = galaxyPredictionsNumbers[badGalaxyPredictions]

starSaveInformation = np.array([badStarIndices, badStarPredictions]).transpose()
galaxySaveInformation = np.array([badGalaxyIndices, badGalaxyPredictions]).transpose()

np.savetxt('data/deconv/starBadPredictionIndexProbability.txt', starSaveInformation)
np.savetxt('data/deconv/galaxyBadPredictionIndexProbability.txt', galaxySaveInformation)

# indices correspond to indices in catalog above
def AnalyzeBadPoints(predictions, indices, magnitudes1, title, xaxistitle, fileName):
    plt.figure()
    color = magnitudes1[indices]
    plt.title(title)
    plt.ylabel(r'$P(Star)$')
    plt.xlabel(xaxistitle)

    qr1x = np.percentile(color, 25)
    qr3x = np.percentile(color, 75)
    IQRx = qr3x - qr1x
    # do it for x
    xticks, xticklabels = plt.xticks()

    if np.min(color) < (qr1x - 1.5 * IQRx):
        xmin = qr1x - 1.5 * IQRx
    else:
        xmin = (3*xticks[0] - xticks[1])/2.
    # shift half a step to the right

    if np.max(color) > (qr3x + 1.5 * IQRx):
        xmax = qr3x + 1.5 * IQRx
    else:
        xmax = (3*xticks[-1] - xticks[-2])/2.

    plt.xlim(xmin, xmax)
    plt.xticks(xticks)

    plt.scatter(color, predictions)
    plt.savefig(fileName)

AnalyzeBadPoints(badStarPredictions, badStarIndices, catalog[:, cold['magR']], \
        'P(Star) for False Negatives in R', 'R (AB)', 'data/PStar/false_negatives_r_star.png')
AnalyzeBadPoints(badStarPredictions, badStarIndices, catalog[:, cold['magI']], \
        'P(Star) for False Negatives in IZ', 'I (AB)', 'data/PStar/false_negatives_i_star.png')
AnalyzeBadPoints(badStarPredictions, badStarIndices, catalog[:, cold['magZ']], \
        'P(Star) for False Negatives in Z', 'Z (AB)', 'data/PStar/false_negatives_z_star.png')

AnalyzeBadPoints(badStarPredictions, badStarIndices, catalog[:, cold['mag_j']], \
        'P(Star) for False Negatives in J', 'J (AB)', 'data/PStar/false_negatives_j_star.png')
AnalyzeBadPoints(badStarPredictions, badStarIndices, catalog[:, cold['mag_h']],  \
        'P(Star) for False Negatives in H', 'H (AB)', 'data/PStar/false_negatives_h_star.png')
AnalyzeBadPoints(badStarPredictions, badStarIndices, catalog[:, cold['mag_k']],  \
        'P(Star) for False Negatives in K', 'H (AB)', 'data/PStar/false_negatives_k_star.png')

AnalyzeBadPoints(badStarPredictions, badStarIndices, catalog[:, cold['mag_36']], \
        r'P(Star) for False Negatives in $3.6\mu{}m$', r'$3.6\mu{}m$ (AB)', \
        'data/PStar/false_negatives_chan1_star.png')
AnalyzeBadPoints(badStarPredictions, badStarIndices, catalog[:, cold['mag_45']], \
        r'P(Star) for False Negatives in $4.5\mu{}m$', r'$4.5\mu{}m$ (AB)', \
        'data/PStar/false_negatives_chan2_star.png')

AnalyzeBadPoints(badGalaxyPredictions, badGalaxyIndices, catalog[:, cold['magR']], \
        'P(Star) for False Positives in R', 'R (AB)', 'data/PStar/false_positives_r_galaxy.png')
AnalyzeBadPoints(badGalaxyPredictions, badGalaxyIndices, catalog[:, cold['magI']], \
        'P(Star) for False Positives in I', 'I (AB)', 'data/PStar/false_positives_i_galaxy.png')
AnalyzeBadPoints(badGalaxyPredictions, badGalaxyIndices, catalog[:, cold['magZ']], \
        'P(Star) for False Positives in Z', 'Z (AB)', 'data/PStar/false_positives_z_galaxy.png')

AnalyzeBadPoints(badGalaxyPredictions, badGalaxyIndices, catalog[:, cold['mag_j']], \
        'P(Star) for False Positives in J', 'J (AB)', 'data/PStar/false_positives_j_galaxy.png')
AnalyzeBadPoints(badGalaxyPredictions, badGalaxyIndices, catalog[:, cold['mag_h']], \
        'P(Star) for False Positives in H', 'H (AB)', 'data/PStar/false_positives_h_galaxy.png')
AnalyzeBadPoints(badGalaxyPredictions, badGalaxyIndices, catalog[:, cold['mag_k']], \
        'P(Star) for False Positives in K', 'K (AB)', 'data/PStar/false_positives_k_galaxy.png')

AnalyzeBadPoints(badGalaxyPredictions, badGalaxyIndices, catalog[:, cold['mag_36']], \
        r'P(Star) for False Positives in $3.6\mu{}m$', r'$3.6\mu{}m$ (AB)', \
        'data/PStar/false_positives_chan1_galaxy.png')
AnalyzeBadPoints(badGalaxyPredictions, badGalaxyIndices, catalog[:, cold['mag_45']], \
        r'P(Star) for False Positives in $4.5\mu{}m$', r'$4.5\mu{}m$ (AB)', \
        'data/PStar/false_positives_chan2_galaxy.png')
