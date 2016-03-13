import itertools
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
def predictStar(clfstar, clfgalaxy, X, Xerr, index):
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

def TryModel(nGaussiansStar, nGaussiansGalaxy): 
    print 'Star Gaussians: {0}'.format(nGaussiansStar)
    print 'Galaxy Gaussians: {0}'.format(nGaussiansGalaxy)

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
    starPredictions = np.array([predictStar(clfstar, clfgalaxy, np.array([XTestStar[i]]), np.array([XErrTestStar[i]]), i) for i in range(starTestNumber)])
    galaxyPredictions = np.array([predictStar(clfstar, clfgalaxy, np.array([XTestGalaxy[i]]), np.array([XErrTestGalaxy[i]]), i) for i in range(galaxyTestNumber)])

    predictions = np.array(starPredictions.tolist() +  galaxyPredictions.tolist())
    results = np.array([1 for i in range(len(starPredictions))] + [0 for i in range(len(galaxyPredictions))])
    report = generateReport(predictions, results)
    return (report['Precision'], clfstar, clfgalaxy)

maxr = None
maxprecision = -1
for r in itertools.product(np.arange(1, 31), np.arange(1, 31)):
    print 'Trying Model Star: {0}, Galaxy: {1}'.format(r[0], r[1])
    precision, clfstar, clfgalaxy = TryModel(r[0], r[1])
    if precision > maxprecision:
        maxprecision = precision
        maxr = r
starPredictions = np.array([predictStar(clfstar, clfgalaxy, np.array([XTestStar[i]]), np.array([XErrTestStar[i]]), i) for i in range(starTestNumber)])

galaxyPredictions = np.array([predictStar(clfstar, clfgalaxy, np.array([XTestGalaxy[i]]), np.array([XErrTestGalaxy[i]]), i) for i in range(galaxyTestNumber)])

predictions = np.array(starPredictions.tolist() +  galaxyPredictions.tolist())
results = np.array([1 for i in range(len(starPredictions))] + [0 for i in range(len(galaxyPredictions))])
report = generateReport(predictions, results)

json.dump(report, open('data/deconv/deconv_results.json', 'w'))

gaussianResults = open('data/deconv/GaussianResults.txt', 'w')
gaussianResults.write('Precision {0}\n'.format(maxprecision))
gaussianResults.write('Gaussians Star {0}\n'.format(maxr[0]))
gaussianResults.write('Gaussians Galaxy {0}\n'.format(maxr[1]))
gaussianResults.close()


#analyze bad classifications
badStarPredictions = starPredictions == 0
badStarIndices = starIndicesTest[badStarPredictions]

badGalaxyPredictions = galaxyPredictions == 1
badGalaxyIndices = galaxyIndicesTest[badGalaxyPredictions]

def DeconvolutionGraphAnalysis(indices, bandColumn, bandTitle, starBoolean):
    magnitudes = catalog[bandColumn][indices]
    plt.figure()
    plt.hist(magnitudes)

    if starBoolean:
        plt.title(r'Distribution of misclassifications for Stars in Band ' + bandTitle)
        plt.xlabel('Magnitude in ' + bandTitle + ' Band')
        plt.ylabel(r'Counts')
        plt.savefig('data/deconv/misclassificationStar' + bandColumn + '.png')
    else:
        plt.title(r'Distribution of misclassifications for Galaxies in Band ' + bandTitle)
        plt.xlabel('Magnitude in ' + bandTitle + ' Band')
        plt.ylabel(r'Counts')
        plt.savefig('data/deconv/misclassificationGalaxy' + bandColumn + '.png')
    
print 'Graphing'
# to be able to use the following function
catalog = pyfits.getdata('MatchHellDeconv2.fits')
DeconvolutionGraphAnalysis(badStarIndices, 'magR', 'R', 1)
DeconvolutionGraphAnalysis(badStarIndices, 'magI', 'I', 1)
DeconvolutionGraphAnalysis(badStarIndices, 'magZ', 'Z', 1)
DeconvolutionGraphAnalysis(badStarIndices, 'mag_j', 'J', 1)
DeconvolutionGraphAnalysis(badStarIndices, 'mag_h', 'H', 1)
DeconvolutionGraphAnalysis(badStarIndices, 'mag_k', 'K', 1)
DeconvolutionGraphAnalysis(badStarIndices, 'mag_36', r'3.6 $\mu{}m$', 1)
DeconvolutionGraphAnalysis(badStarIndices, 'mag_45', r'4.5 $\mu{}m$', 1)

DeconvolutionGraphAnalysis(badGalaxyIndices, 'magR', 'R', 0)
DeconvolutionGraphAnalysis(badGalaxyIndices, 'magI', 'I', 0)
DeconvolutionGraphAnalysis(badGalaxyIndices, 'magZ', 'Z', 0)
DeconvolutionGraphAnalysis(badGalaxyIndices, 'mag_j', 'J', 0)
DeconvolutionGraphAnalysis(badGalaxyIndices, 'mag_h', 'H', 0)
DeconvolutionGraphAnalysis(badGalaxyIndices, 'mag_k', 'K', 0)
DeconvolutionGraphAnalysis(badGalaxyIndices, 'mag_36', r'3.6 $\mu{}m$', 0)
DeconvolutionGraphAnalysis(badGalaxyIndices, 'mag_45', r'4.5 $\mu{}m$', 0)
