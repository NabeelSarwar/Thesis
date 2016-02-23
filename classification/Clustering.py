import json
import sklearn.mixture
import sklearn.svm
import sklearn.linear_model
import matplotlib
import gc
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.markers as markers
plt.ioff()

import numpy as np
import pyfits
#output false positive rate and all that jazz

def generateReport(model, data, results):
    if data.shape[0] != len(results):
        print 'Incorrect size of data and results'
        return {}
    truepositive = 0
    falsepositive = 0
    truenegative = 0
    falsenegative = 0
    predictions = model.predict(data)
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

def tryModel(model, data, ids, results):
    trainingReport = {}
    testReport = {}
    trainingSetSize = int(data.shape[0] * 0.8)
    testSetSize = data.shape[0] - trainingSetSize
    trainingIndices = np.random.choice(data.shape[0], size= trainingSetSize, replace=False)
    allIndices = range(data.shape[0])
    testIndices = filter(lambda x: not x in trainingIndices, allIndices)
    newResults = []
    for e in results:
        if e==2:
            newResults.append(1)
        elif e==1:
            newResults.append(0)
        else:
            print 'Unknown object?'
            print e
            newResults.append(e)
    results = newResults
    results = np.array(results)
    trainingSet = data[trainingIndices, :]
    testSet = data[testIndices, :]
    if testSet.shape[0] != testSetSize:
        print 'Size of test set does not match'
        print testSet.shape[0]
        print testSetSize
    if trainingSet.shape[0] != trainingSetSize:
        print 'Size of training set does not match'
        print trainingSet.shape[0]
        print trainingSetSize

    model.fit(trainingSet, results[trainingIndices])
    trainingReport = generateReport(model, trainingSet, results[trainingIndices])
    testReport = generateReport(model, testSet, results[testIndices])
    return {'training': trainingReport, 'test': testReport}

matchedCat = pyfits.getdata('MatchHell2.fits')
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
goodIndices = magRError < 0.2
goodIndices = sand(goodIndices, magIError < 0.2)
goodIndices = sand(goodIndices, magZError < 0.2)
goodIndices = sand(goodIndices, matchedCat['mag_j_error'] < 0.2)
goodIndices = sand(goodIndices, matchedCat['mag_h_error'] < 0.2)
goodIndices = sand(goodIndices, matchedCat['mag_k_error'] < 0.2)
goodIndices = sand(goodIndices, matchedCat['mag_36error'] < 0.2)
goodIndices = sand(goodIndices, matchedCat['mag_45error'] < 0.2)


# now  we have the proper errors
matchedCat = matchedCat[goodIndices]
print 'Good indices'
print np.any(goodIndices)

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

starindices = np.where(matchedCat['mu_class'] == 2)[0]
galaxyindices = np.where(matchedCat['mu_class'] == 1)[0]

print starindices
print galaxyindices
starGalHash = {}
for e in starindices:
    starGalHash[e] = 'star'

for e in galaxyindices:
    starGalHash[e] = 'galaxy'


results = matchedCat['mu_class']
ids = matchedCat['id']

magnitudeMatrix = np.matrix([magR-magI, magI-magZ, jmags-hmags, hmags-kmags, kmags-chan1mags, chan1mags-chan2mags], dtype=np.dtype('float64')).transpose()


#this dtype somehow encodes the names into the format, so the sklearn thinks that the names are a type and fail to cast
#by names, i mean the column names indicated in the dictionary field names below
#this might be useful later on

#typ = np.dtype({'names': ['magR', 'magI', 'magZ', 'magJ', 'magH', 'magK', 'magChan3', 'magChan2'], 'formats':['float64' for \
#        i in range(magnitudeMatrix.shape[1])]})

#all of this was to be able to restructure the magnitudeMatrix to have floats with column names
#but it was all for naught
typ = np.dtype('float64')
holdMatrix = np.zeros(magnitudeMatrix.shape, dtype=typ)
holdMatrix[:] = magnitudeMatrix[:]
magnitudeMatrix = holdMatrix

numComp = 10
modelGMM = sklearn.mixture.GMM(n_components=numComp)

modelSVM = sklearn.svm.LinearSVC()
logisticModel = sklearn.linear_model.LogisticRegression()


reportSVM = tryModel(modelSVM, magnitudeMatrix, ids, results)
reportLogistic = tryModel(logisticModel, magnitudeMatrix, ids, results)
json.dump(reportSVM, open('data/svm_results.json', 'w'))
json.dump(reportLogistic, open('data/logisitic_regression_results.json', 'w'))
print 'Done'

