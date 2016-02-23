import pyfits
import numpy as np

matchedCat  = pyfits.getdata('MatchHell2.fits')

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
orig_cols = pyfits.ColDefs([pyfits.Column(name = col.name, format=col.format, array=np.array(col.array)) for col in orig_cols])
new_cols = pyfits.ColDefs(cols)
hdu = pyfits.BinTableHDU.from_columns(orig_cols + new_cols)
#remove the file if it exists
try:
    os.remove('MatchHellHSCMags2.fits')
except:
    pass
hdu.writeto('MatchHellHSCMags2.fits')
