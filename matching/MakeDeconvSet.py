import pyfits
import numpy as np
import os

catalog = pyfits.getdata('MatchHellErrorCut2.fits')

columns = ['magR', 'magI', 'magZ', 'mag_j' ,'mag_h', 'mag_k', 'mag_36', 'mag_45']
columns = ['magR', 'magI', 'magZ', 'mag_j' ,'mag_h', 'mag_k', 'mag_36', 'mag_45', 'magRError', 'magIError', 'magZError', 'mag_j_error', \
        'mag_h_error', 'mag_k_error', 'mag_36error', 'mag_45error']


# for an index, says if that catalog has the magnitudes we want
data = np.array([np.logical_or(np.isnan(catalog[col]), np.isinf(catalog[col])) for col in columns]).transpose()
indices = np.array([np.all(np.logical_not(data[i])) for i in range(len(data))])


orig_cols = catalog.columns

orig_cols = pyfits.ColDefs([pyfits.Column(name=col.name, format=col.format, array=np.array(col.array)[indices]) for col in orig_cols])

hdu = pyfits.BinTableHDU.from_columns(orig_cols)

try:
    os.remove('MatchHellDeconv2.fits')
except:
    pass

hdu.writeto('MatchHellDeconv2.fits')
