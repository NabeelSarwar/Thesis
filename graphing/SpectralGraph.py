import pyfits
import numpy as np

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

catalog = pyfits.getdata('MatchHellDeconv2.fits')
# first column is bad index, second is pStar
badStarInfo = np.genfromtxt('data/deconv/starBadPredictionIndexProbability.txt')
magR = catalog['magR']
magRError = catalog['magRError']
magI = catalog['magI']
magIError = catalog['magIError']
magZ = catalog['magZ']
magZError = catalog['magZError']
magY = catalog['magY']
magYError = catalog['magYError']
magJ = catalog['mag_j']
magJError = catalog['mag_j_error']
magH = catalog['mag_h']
magHError = catalog['mag_h_error']
magK = catalog['mag_k']
magKError = catalog['mag_k_error']
mag36 = catalog['mag_36']
mag36Error = catalog['mag_36error']
mag45 = catalog['mag_45']
mag45Error = catalog['mag_45error']


def PlotSpectra(index, filename):
    fig = plt.figure()
    # get current axis
    ax = plt.gca()
    y = [magR[index], magI[index], magZ[index], magY[index], magJ[index], magH[index], magK[index], mag36[index], mag45[index]]
    yerr = [magRError[index], magIError[index], magZError[index], magYError[index], \
            magJError[index], magHError[index], magKError[index], \
            mag36Error[index], mag45Error[index]]

    x = range(len(y))
    plt.errorbar(x, y, yerr = yerr, fmt='.', ecolor='red')
    ax.set_title('Spectra for a Star')
    ax.set_ylabel('Mag (AB)')
    ax.set_xlabel('Band')

    # to get nice ranges
    # this needs to happen after plotting otherwise plot gets messed up

    # do it for x
    xticks, xticklabels = plt.xticks()
    # shift half a step to the left
    # x0 - (x1 - x0) / 2 = (3 * x0 - x1) / 2
    xmin = (3*xticks[0] - xticks[1])/2.
    # shift half a step to the right
    xmax = (3*xticks[-1] - xticks[-2])/2.
    plt.xlim(xmin, xmax)
    plt.xticks(xticks)

    # do it for y
    yticks, yticklabels = plt.yticks()
    # shift half a step to the left
    # x0 - (x1 - x0) / 2 = (3 * x0 - x1) / 2
    ymin = (3*yticks[0] - yticks[1])/2.
    # shaft half a step to the right
    ymax = (3*yticks[-1] - yticks[-2])/2.
    plt.ylim(ymin, ymax)
    plt.yticks(yticks)
    
    ax.set_ylim(ax.get_ylim()[::-1])

    # the first entry is for the zero point
    ax.set_xticklabels(['R', 'I', 'Z', 'Y', 'J', 'H', 'K', 'Chan 1', 'Chan 2'])
    plt.savefig(filename)
    plt.close(fig)


for index in badStarInfo[:, 0]:
    fileName = 'data/spectra/bad_star_{0}.png'.format(int(index))
    PlotSpectra(index, fileName)

# do the same thing for the stars that are actually good, but only do a subset of them
starIndices = np.where(catalog['mu_class']==2)[0].tolist()

for index in badStarInfo[:, 0]:
    if index in starIndices:
        starIndices.remove(index)

starIndices = np.random.choice(np.array(starIndices), 50, replace=False)

for index in starIndices:
    fileName = 'data/spectra/good_star_{0}.png'.format(int(index))
    PlotSpectra(index, fileName)

