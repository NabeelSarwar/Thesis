#needed for ssh graphing
import matplotlib
matplotlib.use('Agg')

from copy import deepcopy
import numpy as np
import pyfits
import matplotlib.pyplot as plt
from astropy.table import Table, Column
from astropy.io import ascii
from sklearn.cluster import KMeans
from astroML.plotting.tools import draw_ellipse
from scipy.stats import gaussian_kde
from sklearn.mixture import GMM

import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom

#needed to graph via ssh
plt.ioff()

def genCatFromTxt(fName, cNames, cNumbers, fTypes=None, cDocs=None, inDegrees=True, **kargs):
    """
    Load data from a `.txt` file into an LSST simple catalog.

    It's assumed that the column names in `cNames` contain the corresponding
    names for the fields of the LSST minimal schema, i.e. `id` `coord.ra` and `coord.dec`.
    """

    if cDocs is None:
        cDocs = ['']*len(cNames)

    if fTypes is None:
        fTypes = ['F']*len(cNames)

    data = np.genfromtxt(fName, names=cNames, usecols=cNumbers, **kargs)

    if inDegrees:
        data['coordra'] = np.radians(data['coordra'])
        data['coorddec'] = np.radians(data['coorddec'])

    schema = afwTable.SimpleTable.makeMinimalSchema()

    for i in range(len(cNames)):
        if cNames[i] not in ['id', 'coord.ra', 'coord.dec']:
            if fTypes[i] != 'String':
                schema.addField(afwTable.Field[fTypes[i]](cNames[i], cDocs[i]))
            else:
                schema.addField(afwTable.Field[fTypes[i]](cNames[i], cDocs[i], 10))

    cat = afwTable.SimpleCatalog(schema)
    n = len(data)
    cat.reserve(n)
    for i in range(n):
        cat.addNew()

    for i, name in enumerate(cNames):
        if fTypes[i] != 'String':
            cat[name][:] = data[name.replace(".", "")]
        else:
            for j, record in enumerate(cat):
                record[name] = data[name.replace(".", "")][j]

    return cat

def genCatFromFits(fName, cNames, fitsNames, fTypes=None, cDocs=None, inDegrees=True, twoDCoord=False):
    """
    Load data from a `.fits` file into an LSST simple catalog.

    It's assumed that the column names in `cNames` contain the corresponding
    names for the fields of the LSST minimal schema, i.e. `id` `coord.ra` and `coord.dec`.
    """

    if cDocs is None:
        cDocs = ['']*len(cNames)

    if fTypes is None:
        fTypes = ['F']*len(cNames)

    hdulist = pyfits.open(fName)
    tbdata = hdulist[1].data

    if inDegrees:
        if twoDCoord:
            # coord must be in position 2
            tbdata[fitsNames[2]][:] = np.radians(tbdata[fitsNames[2]])
        else:
            tbdata[fitsNames[0]][:] = np.radians(tbdata[fitsNames[0]])
            tbdata[fitsNames[1]][:] = np.radians(tbdata[fitsNames[1]])


    schema = afwTable.SimpleTable.makeMinimalSchema()

    for i in range(len(cNames)):
        if cNames[i] not in ['id', 'coord.ra', 'coord.dec', 'coord'] and cNames[i] not in schema:
            schema.addField(afwTable.Field[fTypes[i]](cNames[i], cDocs[i]))

    cat = afwTable.SimpleCatalog(schema)
    n = len(tbdata)
    cat.reserve(n)
    for i in range(n):
        cat.addNew()

    for i in range(len(cNames)):
        if twoDCoord and not cNames[i] in ['coord', 'coord.ra', 'coord.dec']: # #and (cNames[i] in ['coord.ra', 'coord.dec'])):
            cat[cNames[i]][:] = tbdata[fitsNames[i]]
        elif not twoDCoord:
            cat[cNames[i]][:] = tbdata[fitsNames[i]]
    if twoDCoord:
        cat['coord.ra'][:] = tbdata[fitsNames[1]][:,0]
        cat['coord.dec'][:] = tbdata[fitsNames[1]][:,1]

    return cat

def genCatFromIpac(fName, cNames, ipacNames, fTypes=None, cDocs=None, inDegrees=True, twoDCoord=False):
    """
    Load data from a table in ipac format into an LSST simple catalog.

    It's assumed that the column names in `cNames` contain the corresponding
    names for the fields of the LSST minimal schema, i.e. `id` `coord.ra` and `coord.dec`.
    """

    if cDocs is None:
        cDocs = ['']*len(cNames)

    if fTypes is None:
        fTypes = ['F']*len(cNames)

    table = ascii.read(fName)
    #table = Table.read(fName, format='ipac')

    if inDegrees:
        if twoDCoord:
            table[ipacNames[1]].data.data[:] = np.radians(table[ipacNames[1]].data.data)
        else:
            table[ipacNames[0]] = np.radians(table[ipacNames[0]])
            table[ipacNames[1]] = np.radians(table[ipacNames[1]])

    schema = afwTable.SimpleTable.makeMinimalSchema()

    for i in range(len(cNames)):
        if cNames[i] not in ['id', 'coord.ra', 'coord.dec', 'coord']:
            schema.addField(afwTable.Field[fTypes[i]](cNames[i], cDocs[i]))

    cat = afwTable.SimpleCatalog(schema)
    n = len(table)
    cat.reserve(n)
    for i in range(n):
        cat.addNew()

    for i in range(len(cNames)):
        if not (twoDCoord and (cNames[i] in ['coord.ra', 'coord.dec'])):
            try:
                cat[cNames[i]][:] = table[ipacNames[i]].data
            except ValueError:
                import ipdb; ipdb.set_trace()
    if twoDCoord:
        cat['coord.ra'][:] = tbdata[ipacNames[1]].data.data[:,0]
        cat['coord.dec'][:] = tbdata[ipacNames[1]].data.data[:,1]

    return cat

def matchCats(cat1, cat2, matchRadius=1*afwGeom.arcseconds, includeMismatches=True, multiMeas=False):
    """
    Match to catalogs and return a catalog with the fields of the two catalogs
    """

    mc = afwTable.MatchControl()
    mc.includeMismatches = includeMismatches
    mc.findOnlyClosest = True

    matched = afwTable.matchRaDec(cat1, cat2, matchRadius, mc)

    bestMatches = {}
    if includeMismatches:
        noMatch = []
    for m1, m2, d in matched:
        if m2 is None:
            noMatch.append(m1)
        else:
            if not multiMeas:
                id = m2.getId()
                if id not in bestMatches:
                    bestMatches[id] = (m1, m2, d)
                else:
                    if d < bestMatches[id][2]:
                        bestMatches[id] = (m1, m2, d)
            else:
                id = m1.getId()
                bestMatches[id] = (m1, m2, d)

    if includeMismatches:
        print "{0} objects from {1} in the first catalog had no match in the second catalog.".format(len(noMatch), len(cat1))
        print "{0} objects from the first catalog with a match in the second catalog were not the closest match.".format(len(matched) - len(noMatch) - len(bestMatches))

    nMatches = len(bestMatches)
    print "I found {0} matches".format(nMatches)

    schema1 = cat1.getSchema(); schema2 = cat2.getSchema()
    names1 = cat1.schema.getNames(); names2 = cat2.schema.getNames()

    schema = afwTable.SimpleTable.makeMinimalSchema()

    catKeys = []; cat1Keys = []; cat2Keys = []
    for name in names1:
        cat1Keys.append(schema1.find(name).getKey())
        if name not in ['id', 'coord']:
            catKeys.append(schema.addField(schema1.find(name).getField()))
        else:
            catKeys.append(schema.find(name).getKey())
    for name in names2:
        cat2Keys.append(schema2.find(name).getKey())
        if name not in schema1.getNames():
            catKeys.append(schema.addField(schema2.find(name).getField()))
        elif name+".2" not in schema1.getNames():
            catKeys.append(schema.addField(schema2.find(name).getField().copyRenamed(name+".2")))
        else:
            text = 3
            if not name +'.'+str(text) in schema.getNames():
                key = schema.addField(schema2.find(name).getField().copyRenamed(name+"." +str(text)))
                catKeys.append(key)
            else:
                while name+'.'+str(text) in schema.getNames():
                    print 'insane looping 2:' + name+'.'+str(text)
                    text = text + 1
                    if not name+'.'+str(text) in schema.getNames():
                        catKeys.append(schema.addField(schema2.find(name).getField().copyRenamed(name+"."+str(text))))
                        # need to break or otherwise introduce infinite loop because
                        # the thing will be added and then we try to add next iteration
                        break

    print 'Done matching'
    print 'Now merging'
    cat = afwTable.SimpleCatalog(schema)
    cat.reserve(nMatches)

    for id in bestMatches:
        m1, m2, d = bestMatches[id]
        record = cat.addNew()
        for i in range(len(cat1Keys)):
            record.set(catKeys[i], m1.get(cat1Keys[i]))
        for i in range(len(cat1Keys), len(catKeys)):
            record.set(catKeys[i], m2.get(cat2Keys[i-len(cat1Keys)]))

    return cat

print 'Loading HSC'
hscdata = pyfits.getdata('matchDeepCoaddMeas-137520151126CosmosGRIZY.fits')
columns = np.array([w.name for w in hscdata.columns])
shape = np.array([len(hscdata[w].shape) for w in columns])
columns = columns[np.logical_or(shape == 1, columns=='coord')].tolist()
hschst = genCatFromFits('matchDeepCoaddMeas-137520151126CosmosGRIZY.fits', columns, columns, inDegrees=False, twoDCoord=True)

print 'Loading Vista'
uvistah = genCatFromFits('UVISTA_Ks_01_09_13_allpaw_skysub_015_dr2_v2_H.cat.fits', ['coord.ra', 'coord.dec', 'mag_h', 'mag_h_error'], ['ALPHA_J2000', 'DELTA_J2000', 'MAG_AUTO', 'MAGERR_AUTO'])
uvistaj = genCatFromFits('UVISTA_Ks_01_09_13_allpaw_skysub_015_dr2_v2_J.cat.fits', ['coord.ra', 'coord.dec', 'mag_j', 'mag_j_error'], ['ALPHA_J2000', 'DELTA_J2000', 'MAG_AUTO', 'MAGERR_AUTO'])
uvistak = genCatFromFits('UVISTA_Ks_01_09_13_allpaw_skysub_015_dr2_v2_Ks.cat.fits', ['coord.ra', 'coord.dec', 'mag_k', 'mag_k_error'], ['ALPHA_J2000', 'DELTA_J2000', 'MAG_AUTO', 'MAGERR_AUTO'])
uvistay = genCatFromFits('UVISTA_Ks_01_09_13_allpaw_skysub_015_dr2_v2_Y.cat.fits', ['coord.ra', 'coord.dec', 'mag_y', 'mag_y_error'], ['ALPHA_J2000', 'DELTA_J2000', 'MAG_AUTO', 'MAGERR_AUTO'])

print 'Ipac stuff'
#ipac information
# need to change name and unit of column and do this for all different of the filters
spitzerData = ascii.read('SpitzerData.tbl')

names = [[str(fil), 'mag_' + str(fil) + '_error'] for fil in [36, 45, 58, 80]]
names = sum(names, [])
# get proper indices
spitzerData = spitzerData[np.logical_and.reduce((spitzerData['flux_c1_2'] > 0, spitzerData['flux_c2_2'] > 0 , spitzerData['flux_c3_2'] > 0, spitzerData['flux_c4_2'] > 0 ))]

# getting the magnitudes using the spitzer data
mag36 = Column(data=-2.5 * np.log10(spitzerData['flux_c1_2'].data/0.765) + 23.9 - 2.788, name ='mag_36', unit='AB mag')
mag45 = Column(data=-2.5 * np.log10(spitzerData['flux_c2_2'].data/0.740) + 23.9 -3.255, name ='mag_45', unit='AB mag')
mag58 = Column(data=-2.5 * np.log10(spitzerData['flux_c3_2'].data/0.625) + 23.9 -3.743, name = 'mag_58', unit = 'AB mag')
mag80 = Column(data=-2.5 * np.log10(spitzerData['flux_c4_2'].data/0.580) + 23.9 -4.372, name = 'mag_80', unit = 'AB mag')

#include magnitude errors
mag36error = Column(data = np.abs(1.08574/spitzerData['flux_c1_2'].data*spitzerData['err_c1_2'].data), name='mag_36error', unit='AB mag')
mag45error = Column(data = np.abs(1.08574/spitzerData['flux_c2_2'].data * spitzerData['err_c2_2'].data), name='mag_45error', unit='AB mag')
mag58error = Column(data = np.abs(1.08574/ spitzerData['flux_c3_2'].data*spitzerData['err_c3_2'].data), name='mag_58error', unit='AB mag')
mag80error = Column(data = np.abs(1.08574/spitzerData['flux_c4_2'].data *spitzerData['err_c4_2'].data), name='mag_80error', unit='AB mag')
spitzerData.add_columns([mag36, mag45, mag58, mag80, mag36error, mag45error, mag58error, mag80error])
ascii.write(output='SpitzerWithMag.tbl', table=spitzerData, format='ipac')

spitzerCatalog = genCatFromIpac('SpitzerWithMag.tbl', ['coord.ra', 'coord.dec', 'mag_36', 'mag_45', 'mag_58', 'mag_80', 'mag_36error', 'mag_45error', 'mag_58error', 'mag_80error'], \
['ra', 'dec', 'mag_36', 'mag_45', 'mag_58', 'mag_80', 'mag_36error', 'mag_45error', 'mag_58error', 'mag_80error'])


matchedCat = matchCats(uvistah, uvistaj)
matchedCat = matchCats(matchedCat, uvistak)
matchedCat = matchCats(matchedCat, uvistay)

uvistahsc = matchCats(matchedCat, hschst)
# hsc and vista match
uvistahsc.writeFits('UltraVistaHSCHST.fits')

# this is the overall
matchedCat = matchCats(matchedCat, spitzerCatalog)

# hsc and spitzer match
spitzerHSC = matchCats(hschst, spitzerCatalog)
spitzerHSC.writeFits('SpitzerHSCHST.fits')

# load in jose's work
# thus finish overall matching
print 'I got to here'
matchedCat.writeFits('MatchHell2.fits')
