import numpy as np
import astropy.io.fits as fits
from astropy.table import Table

sr_masks = ['masks/mask_hi_25.fits', 'masks/mask_hi_18.fits', 'masks/mask_hi_08.fits']
outnames = ['hi_fir_gt18M31data.csv','hi_fir_bet818M31data.csv','hi_fir_lt8M31data.csv']

def doall(sfr='sfr_fir_reg_hi.fits',gas='hi_per_mass_sun_pc2.fits',masklist=sr_masks):
    sfrdat = fits.getdata(sfr)
    gasdat = fits.getdata(gas)
    sfr_unc = fits.getdata('err_'+sfr)
    gas_unc = fits.getdata('err_'+gas)

    for i,mask in enumerate(masklist):
        position_mask = fits.getdata(mask)
        data_mask = np.logical_and(gasdat>0,np.logical_not(np.isnan[sfrdat]))
        fit_mask = np.logical_and(~np.isnan(position_mask), data_mask) 
        gas = np.log10(gasdat[fit_mask])
        sfr = np.log10(sfrdat[fit_mask])
        gas_unc = np.abs(gas_unc[fit_mask]/(np.log(10)*gasdat[fit_mask]))
        sfr_unc = np.abs(sfr_unc[fit_mask]/(np.log(10)*sfrdat[fit_mask]))
        index = np.ones(len(gas)).astype(int)*31
        fittab = Table([gas, gas_unc, sfr, sfr_unc])
        fittab.write(outnames[i], format='ascii.csv')
    return
