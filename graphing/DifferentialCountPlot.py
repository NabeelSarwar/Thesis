import pyfits
import matplotlib.pyplot as plt
import numpy as np
import matplotlib

matplotlib.use('Agg')
catalog = pyfits.getdata('MatchHellErrorCut2.fits')

def DifferentialCountPlot(starBoolean, catalog, bandColumn, bandTitle, bins, fileName):
    if starBoolean:
        entries = catalog[catalog['mu_class']==2]
    else:
        entries = catalog[catalog['mu_class']==1]

    binCounts = np.zeros(len(bins))
    bins2 = np.append(bins, [np.infty])

    magnitudes = np.array(entries[bandColumn])
    magnitudes = magnitudes[np.logical_not(np.isnan(magnitudes))]
    print len(magnitudes)
    for i in range(len(bins2)-1):
        binCounts[i] = np.sum(np.logical_and(magnitudes >= bins2[i], magnitudes < bins2[i+1]))

    fig = plt.figure()
    if starBoolean:
        extra = 'For Stars'
    else:
        extra = 'For Galaxies'

    plt.title(r'Logarithm of Approximate $\frac{dN}{dm}$ in ' + bandTitle + ' ' + extra)
    plt.xlabel('Magnitude in ' + bandTitle + ' Band (AB)')
    plt.ylabel(r'$\log_{10}(\Delta{}N)$')
    binCounts = np.log10(binCounts)
    width = 1.0
    plt.xticks(bins+width, bins)
    plt.bar(bins, binCounts, width)
    plt.savefig(fileName)

bins = np.arange(30)

DifferentialCountPlot(1, catalog, 'magR', 'R', bins, 'data/counts/starsRcounts.png')
DifferentialCountPlot(1, catalog, 'magI', 'I', bins, 'data/counts/starsIcounts.png')
DifferentialCountPlot(1, catalog, 'magZ', 'Z', bins, 'data/counts/starsZcounts.png')
DifferentialCountPlot(1, catalog, 'magY', 'Y', bins, 'data/counts/starsYcounts.png')
DifferentialCountPlot(1, catalog, 'mag_j', 'J', bins, 'data/counts/starsJcounts.png')
DifferentialCountPlot(1, catalog, 'mag_h', 'H', bins, 'data/counts/starsHcounts.png')
DifferentialCountPlot(1, catalog, 'mag_k', 'K', bins, 'data/counts/starsKcounts.png')
DifferentialCountPlot(1, catalog, 'mag_36', r'3.6 $\mu$m', bins, 'data/counts/stars36counts.png')
DifferentialCountPlot(1, catalog, 'mag_45', r'4.5 $\mu$m', bins, 'data/counts/stars45counts.png')

DifferentialCountPlot(0, catalog, 'magR', 'R', bins, 'data/counts/galaxiesRcounts.png')
DifferentialCountPlot(0, catalog, 'magI', 'I', bins, 'data/counts/galaxiesIcounts.png')
DifferentialCountPlot(0, catalog, 'magZ', 'Z', bins, 'data/counts/galaxiesZcounts.png')
DifferentialCountPlot(0, catalog, 'magY', 'Y', bins, 'data/counts/galaxiesYcounts.png')
DifferentialCountPlot(0, catalog, 'mag_j', 'J', bins, 'data/counts/galaxiesJcounts.png')
DifferentialCountPlot(0, catalog, 'mag_h', 'H', bins, 'data/counts/galaxiesHcounts.png')
DifferentialCountPlot(0, catalog, 'mag_k', 'K', bins, 'data/counts/galaxiesKcounts.png')
DifferentialCountPlot(0, catalog, 'mag_36', r'3.6 $\mu$m', bins, 'data/counts/galaxies36counts.png')
DifferentialCountPlot(0, catalog, 'mag_45', r'4.5 $\mu$m', bins, 'data/counts/galaxies45counts.png')
