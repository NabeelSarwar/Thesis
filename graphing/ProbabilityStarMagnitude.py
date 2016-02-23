import pyfits
import numpy as np
import extreme_deconvolution as xd
from sklearn.mixture import GMM
from astroML.density_estimation import XDGMM
import json

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

catalog = pyfits.getdata('MatchHellDeconv2.fits')
columns = ['magR', 'magI', 'magZ', 'mag_j' ,'mag_h', 'mag_k', 'mag_36', 'mag_45', 'magRError', 'magIError', 'magZError', 'mag_j_error', \
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
    colors.append(data[cold['mag_j']] - data[cold['mag_h']])
    colors.append(data[cold['mag_h']] - data[cold['mag_k']])
    colors.append(data[cold['mag_36']] - data[cold['mag_45']])
    colors = np.array(colors)

    if np.any(np.isnan(colors)):
        print index
    return colors 

def  makeNoiseMatrix(index):
    noise = np.zeros(25).reshape(5, 5)

    data = catalog[index]

    noise[0, 0] = data[cold['magRError']]**2 + data[cold['magIError']]**2
    noise[1, 1] = data[cold['magIError']]**2 + data[cold['magZError']]**2
    noise[2, 2] = data[cold['mag_j_error']]**2 + data[cold['mag_h_error']]**2
    noise[3, 3] = data[cold['mag_h_error']]**2 + data[cold['mag_k_error']]**2
    noise[4, 4] = data[cold['mag_36error']]**2 + data[cold['mag_45error']]**2
    noise[0, 1] = -data[cold['magIError']]**2
    noise[2, 3] = -data[cold['mag_h_error']]**2

    noise[1, 0] = noise[0, 1]
    noise[3, 2] = noise[2, 3]

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

nGaussiansStar = 18
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
json.dump(report, open('deconv_results.json', 'w'))


#analyze bad classifications
badStarPredictions = np.where(starPredictions == 0)[0]
badStarIndices = starIndicesTest[badStarPredictions]

badGalaxyPredictions = np.where(galaxyPredictions == 1)[0]
badGalaxyIndices = galaxyIndicesTest[badGalaxyPredictions]

starPredictionsNumbers = np.array([predictStarDouble(np.array([XTestStar[i]]), np.array([XErrTestStar[i]]), starIndicesTest[i]) for i in range(starTestNumber)])

galaxyPredictionsNumbers = np.array([predictStarDouble(np.array([XTestGalaxy[i]]), np.array([XErrTestGalaxy[i]]), galaxyIndicesTest[i]) for i in range(galaxyTestNumber)])

badStarPredictions = starPredictionsNumbers[badStarPredictions]
badGalaxyPredictions = galaxyPredictionsNumbers[badGalaxyPredictions]


# indices correspond to indices in catalog above
def AnalyzeBadPoints(predictions, indices, magnitudes1, magnitudes2, title, xaxistitle, fileName):
    plt.figure()
    color = magnitudes1[indices]-magnitudes2[indices]
    plt.title(title)
    plt.ylabel(r'$P(Star)$')
    plt.xlabel(xaxistitle)

    plt.scatter(color, predictions)
    plt.savefig(fileName)


AnalyzeBadPoints(badStarPredictions, badStarIndices, catalog[cold['magR']], catalog[cold['magI']], \
        'P(Star) for False Negatives in R-I', 'R-I (AB)', 'data/PStar/false_negatives_ri.png')
