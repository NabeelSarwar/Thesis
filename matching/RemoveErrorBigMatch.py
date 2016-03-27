#After making a huge catalog in matchingCode.py, I remove all samples that have any part with magnitude errors larger than 5sigma
#this also adds the columns magG, magR, magI, magZ, magY and their errors
# I had to run parts of this line by line in ipython to get it working, so near the end there could be errors
import pyfits
import numpy as np
import os

matchedCat = pyfits.getdata('MatchHell2.fits')
#shortcut and
sand = np.logical_and

# not the real good indices, but these just get non zero flxues out
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

#use this matchedCat for the proper compilation of matches
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

# add these columns to the table
cols = [] 
cols.append(
    pyfits.Column(name='magR', format='D', array= magR)
    )
cols.append(
    pyfits.Column(name='magI', format='D', array= magI)
    )
cols.append(
    pyfits.Column(name='magG', format='D', array= magG)
    )
cols.append(
    pyfits.Column(name='magZ', format='D', array= magZ)
    )
cols.append(
    pyfits.Column(name='magY', format='D', array= magY)
    )
cols.append(
    pyfits.Column(name='magGError', format='D', array= magGError)
    )
cols.append(
    pyfits.Column(name='magRError', format='D', array= magRError)
    )
cols.append(
    pyfits.Column(name='magIError', format='D', array= magIError)
    )
cols.append(
    pyfits.Column(name='magZError', format='D', array= magZError)
    )
cols.append(
    pyfits.Column(name='magYError', format='D', array= magYError)
    )

#does not change the underlying data
orig_cols = matchedCat.columns
orig_cols = pyfits.ColDefs([pyfits.Column(name = col.name, format=col.format, array=np.array(col.array)[goodIndices]) for col in orig_cols])
new_cols = pyfits.ColDefs(cols)
hdu = pyfits.BinTableHDU.from_columns(orig_cols + new_cols)
os.remove('MatchHellErrorCut2.fits')
hdu.writeto('MatchHellErrorCut2.fits')

#dirty way to not have to deal with the data
matchedCat = pyfits.getdata('MatchHellErrorCut2.fits')

# establish the indices with proper colors, and then get this

#possible errors ahead
goodIndices = magRError < 0.2
goodIndices = sand(goodIndices, magIError < 0.2)
goodIndices = sand(goodIndices, magZError < 0.2)
goodIndices = sand(goodIndices, magGError < 0.2)
goodIndices = sand(goodIndices, magYError < 0.2)
goodIndices = sand(goodIndices, matchedCat['mag_j_error'] < 0.2)
goodIndices = sand(goodIndices, matchedCat['mag_h_error'] < 0.2)
goodIndices = sand(goodIndices, matchedCat['mag_k_error'] < 0.2)
goodIndices = sand(goodIndices, matchedCat['mag_36error'] < 0.2)
goodIndices = sand(goodIndices, matchedCat['mag_45error'] < 0.2)
#possible errors end here

# now  we have the proper errors
orig_cols = matchedCat.columns
orig_cols = pyfits.ColDefs([pyfits.Column(name = col.name, format=col.format, array=np.array(col.array)[goodIndices]) for col in orig_cols])
hdu = pyfits.BinTableHDU.from_columns(orig_cols)
try:
    os.remove('MatchHellErrorCut2.fits')
except:
    pass
hdu.writeto('MatchHellErrorCut2.fits')
